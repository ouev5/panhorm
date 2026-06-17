#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recompute RAG 200 evaluation with an evidence-anchored deterministic rubric.
Uses saved per-question RAG/direct answers from evaluation/rag_200_eval/results.
Does not rewrite original results; writes outputs to evaluation/rag_200_eval_v2.
"""
import json, re, csv, math, statistics
from pathlib import Path
from datetime import datetime

BASE = Path('/www/wwwroot/default/animal_hormone')
SRC = BASE / 'evaluation' / 'rag_200_eval'
OUT = BASE / 'evaluation' / 'rag_200_eval_v2'
FIG = OUT / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

def now(): return datetime.now().isoformat(timespec='seconds')
def norm(s): return re.sub(r'[^a-z0-9]+',' ',str(s or '').lower()).strip()
def words(s): return [w for w in norm(s).split() if len(w) >= 2]
def clamp(x): return max(0.0, min(1.0, float(x)))

def expected_terms(expected):
    raw = str(expected or '').strip()
    if not raw: return []
    parts = re.split(r'[;,|/]+|\band\b', raw, flags=re.I)
    terms=[]
    for p in parts:
        p=p.strip()
        if not p: continue
        # keep compact gene-like symbol or meaningful phrase
        if len(p) <= 80:
            terms.append(p)
    # unique preserving order
    out=[]; seen=set()
    for t in terms:
        k=norm(t)
        if k and k not in seen:
            out.append(t); seen.add(k)
    return out[:8]

def contains_term(answer, term):
    a = norm(answer); t = norm(term)
    if not t: return False
    # exact normalized phrase
    if re.search(r'(^| )' + re.escape(t) + r'( |$)', a): return True
    # gene/protein symbols: exact original case-insensitive word boundary
    if re.fullmatch(r'[A-Za-z0-9_.-]{2,20}', str(term).strip()):
        return re.search(r'\b' + re.escape(str(term).strip()) + r'\b', str(answer or ''), re.I) is not None
    # phrase fallback: most tokens present
    tw = words(term)
    if not tw: return False
    present=sum(1 for w in tw if re.search(r'(^| )'+re.escape(w)+r'( |$)', a))
    return present >= max(1, math.ceil(len(tw)*0.65))

def term_coverage(answer, terms):
    if not terms: return 0.5
    return sum(1 for t in terms if contains_term(answer, t)) / len(terms)

def question_focus_coverage(question, answer):
    # Remove generic words, measure whether answer mentions query entities.
    stop=set('which what is are in the database gene genes hormone associated with for summarize receptor coding information recorded relation disease phenotype and or of to a an by from'.split())
    qwords=[w for w in words(question) if w not in stop and len(w)>2]
    if not qwords: return 0.75 if answer else 0
    a=norm(answer)
    present=sum(1 for w in qwords[:12] if re.search(r'(^| )'+re.escape(w)+r'( |$)', a))
    return clamp(present / min(len(qwords),12))

def answer_quality(answer):
    if not answer: return 0
    low=str(answer).lower()
    if 'error ' in low[:80] or 'request timeout' in low[:120]: return 0
    length_score=clamp(len(answer)/350)
    structure=0.15 if any(x in answer for x in ['-', '•', ':', '\n']) else 0
    caveat=0.10 if any(x in low for x in ['database', 'record', 'provided context', 'based on']) else 0
    return clamp(0.75*length_score + structure + caveat)

def traceability(answer):
    score=0.0
    pats=[r'PMID[:\s]*\d+', r'UniProt|NCBI|PubChem|DOID|DO:', r'gene|receptor|hormone|species|organism', r'database|record|evidence|source|provided context', r'\b[A-Z][A-Za-z0-9_.-]{2,12}\b']
    weights=[0.25,0.15,0.20,0.25,0.15]
    for pat,w in zip(pats,weights):
        if re.search(pat, answer or '', re.I): score+=w
    return clamp(score)

def score_record(d):
    q=d.get('question',''); exp=d.get('expected',''); rag=d.get('rag_answer',''); direct=d.get('direct_answer','')
    terms=expected_terms(exp)
    cov=term_coverage(rag, terms)
    direct_cov=term_coverage(direct, terms)
    focus=question_focus_coverage(q, rag)
    qual=answer_quality(rag)
    tr=traceability(rag)
    # Evidence-anchored metrics. Coverage dominates accuracy; focus dominates relevance.
    relevance=clamp(0.60*focus + 0.25*qual + 0.15*(1 if rag else 0))
    accuracy=clamp(0.82*cov + 0.10*focus + 0.08*tr)
    completeness=clamp(0.55*cov + 0.25*qual + 0.20*tr)
    overall=clamp(0.4*relevance + 0.3*completeness + 0.3*accuracy)
    f1=clamp(2*accuracy*completeness/(accuracy+completeness)) if (accuracy+completeness)>0 else 0
    # Hallucination proxy: missing expected evidence and low question focus increase hallucination risk.
    rag_hall=clamp(0.75*(1-cov) + 0.25*(1-focus))
    direct_hall=clamp(0.80*(1-direct_cov) + 0.20*(1-question_focus_coverage(q,direct)))
    return {
        'expected_terms':terms,
        'expected_coverage':cov,
        'direct_expected_coverage':direct_cov,
        'relevance':relevance,
        'completeness':completeness,
        'accuracy':accuracy,
        'overall':overall,
        'traceability':tr,
        'f1':f1,
        'rag_hallucination':rag_hall,
        'direct_hallucination':direct_hall,
        'passed': overall >= 0.75,
        'focus_coverage':focus,
        'answer_quality':qual,
    }

def p_value_two_rates(a_success, a_n, b_success, b_n):
    # Two-proportion z-test normal approximation.
    if a_n == 0 or b_n == 0: return None
    p1=a_success/a_n; p2=b_success/b_n; p=(a_success+b_success)/(a_n+b_n)
    se=math.sqrt(p*(1-p)*(1/a_n+1/b_n)) if 0 < p < 1 else 0
    if se == 0: return 1.0
    z=(p1-p2)/se
    # two-sided p using erfc
    return math.erfc(abs(z)/math.sqrt(2))

def aggregate(scored):
    done=scored
    def mean(k): return sum(x['scores'][k] for x in done)/len(done) if done else 0
    bytype={}
    for t in sorted(set(x['type'] for x in done)):
        xs=[x for x in done if x['type']==t]
        bytype[t]={'n':len(xs),'f1':sum(x['scores']['f1'] for x in xs)/len(xs),'overall':sum(x['scores']['overall'] for x in xs)/len(xs)}
    rag_hall_flags=sum(1 for x in done if x['scores']['rag_hallucination'] >= 0.5)
    direct_hall_flags=sum(1 for x in done if x['scores']['direct_hallucination'] >= 0.5)
    p=p_value_two_rates(rag_hall_flags, len(done), direct_hall_flags, len(done))
    return {
        'generated_at':now(),
        'rubric':'v2 evidence-anchored deterministic recomputation from saved answers',
        'n_total':200,'n_done':len(done),
        'overall_score':mean('overall'),'relevance':mean('relevance'),'completeness':mean('completeness'),'accuracy':mean('accuracy'),
        'pass_threshold':0.75,'pass_count':sum(1 for x in done if x['scores']['passed']),'pass_rate':sum(1 for x in done if x['scores']['passed'])/len(done),
        'rag_hallucination_rate':mean('rag_hallucination'),'direct_hallucination_rate':mean('direct_hallucination'),
        'rag_hallucination_flag_rate':rag_hall_flags/len(done),'direct_hallucination_flag_rate':direct_hall_flags/len(done),'hallucination_flag_p_value':p,
        'traceability':mean('traceability'),'by_type':bytype,
    }

def make_figures(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax=plt.subplots(figsize=(7,5), dpi=300)
    labels=['RAG hallucination\n(mean risk)','Direct DeepSeek\n(mean risk)','Traceability']
    vals=[summary['rag_hallucination_rate'],summary['direct_hallucination_rate'],summary['traceability']]
    bars=ax.bar(labels, vals, color=['#2ca25f','#de2d26','#3182bd'])
    ax.set_ylim(0,1); ax.set_ylabel('Score / Rate'); ax.set_title('Figure 6C. Hallucination reduction and traceability')
    ax.grid(axis='y', alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}', ha='center', fontsize=9)
    plt.tight_layout(); fig.savefig(FIG/'figure_6C_hallucination_traceability.png'); plt.close(fig)
    types=list(summary['by_type'].keys()); vals=[summary['by_type'][t]['f1'] for t in types]
    fig, ax=plt.subplots(figsize=(7,5), dpi=300)
    bars=ax.bar(types, vals, color=['#756bb1','#31a354','#fd8d3c','#6baed6'][:len(types)])
    ax.set_ylim(0,1); ax.set_ylabel('F1 score'); ax.set_title('Figure 6D. Performance by question type'); ax.grid(axis='y', alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}', ha='center', fontsize=9)
    plt.xticks(rotation=20, ha='right'); plt.tight_layout(); fig.savefig(FIG/'figure_6D_f1_by_question_type.png'); plt.close(fig)
    fig, ax=plt.subplots(figsize=(8,5), dpi=300)
    labels=['Overall','Relevance','Completeness','Accuracy','Pass rate']
    vals=[summary['overall_score'],summary['relevance'],summary['completeness'],summary['accuracy'],summary['pass_rate']]
    bars=ax.bar(labels, vals, color='#4c78a8')
    ax.set_ylim(0,1); ax.set_ylabel('Score'); ax.set_title('RAG 200-question evaluation (recomputed rubric)'); ax.grid(axis='y', alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}', ha='center', fontsize=9)
    plt.tight_layout(); fig.savefig(FIG/'rag_200_overall_metrics.png'); plt.close(fig)

records=[]
for p in sorted((SRC/'results').glob('result_*.json')):
    d=json.loads(p.read_text(encoding='utf-8'))
    if d.get('status')!='done': continue
    s=score_record(d)
    out={k:d.get(k) for k in ['id','type','question','expected','rag_answer','direct_answer','rag_context_type','rag_seconds','direct_seconds']}
    out['scores']=s
    records.append(out)
    (OUT/'results').mkdir(exist_ok=True)
    (OUT/'results'/p.name).write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
summary=aggregate(records)
(OUT/'rag_200_summary_v2.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
with open(OUT/'rag_200_results_v2.csv','w',newline='',encoding='utf-8') as f:
    w=csv.writer(f)
    w.writerow(['id','type','question','expected','overall','relevance','completeness','accuracy','f1','traceability','rag_hallucination','direct_hallucination','passed','expected_coverage','direct_expected_coverage','rag_context_type'])
    for r in records:
        s=r['scores']; w.writerow([r['id'],r['type'],r['question'],r['expected'],s['overall'],s['relevance'],s['completeness'],s['accuracy'],s['f1'],s['traceability'],s['rag_hallucination'],s['direct_hallucination'],s['passed'],s['expected_coverage'],s['direct_expected_coverage'],r.get('rag_context_type')])
make_figures(summary)
print(json.dumps(summary, indent=2, ensure_ascii=False))
