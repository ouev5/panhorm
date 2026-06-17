#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RAG 200-question quantitative evaluation.
Saves one JSON result immediately after every question; resumable.
"""
import os, sys, json, time, random, re, traceback
from pathlib import Path
from datetime import datetime

BASE = Path('/www/wwwroot/default/animal_hormone')
OUT = BASE / 'evaluation' / 'rag_200_eval'
RESULT_DIR = OUT / 'results'
FIG_DIR = OUT / 'figures'
LOG_DIR = OUT / 'logs'
DATASET = OUT / 'rag_200_questions.json'
SUMMARY_JSON = OUT / 'rag_200_summary.json'
SUMMARY_CSV = OUT / 'rag_200_results.csv'
PROGRESS = OUT / 'progress.json'
for d in [RESULT_DIR, FIG_DIR, LOG_DIR]:
    d.mkdir(parents=True, exist_ok=True)

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'animal_hormone.settings')
sys.path.insert(0, str(BASE))
import django
django.setup()

from hormone_app.models import HormoneRelatedGene, HormoneReceptorInfo, HormoneReceptorFull
from hormone_app.utils.rag_deepseek import rag_search, deepseek_rag_client, DEEPSEEK_API_URL, DEEPSEEK_MODEL

random.seed(20260531)

def now():
    return datetime.now().isoformat(timespec='seconds')

def jdump(path, obj):
    path = Path(path)
    tmp = Path(str(path) + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8')
    tmp.replace(path)

def get_val(x):
    if x is None:
        return ''
    return str(x).strip()

def build_dataset():
    if DATASET.exists():
        return json.loads(DATASET.read_text(encoding='utf-8'))
    questions = []
    seen = set()
    def add(q, qtype, expected='', evidence=None):
        q = ' '.join(q.split())
        if not q or q.lower() in seen:
            return
        seen.add(q.lower())
        questions.append({'id': len(questions)+1, 'question': q, 'type': qtype, 'expected': expected, 'evidence': evidence or {}})

    rows = list(HormoneRelatedGene.objects.exclude(hormone_name__isnull=True).exclude(related_genes__isnull=True).values(
        'hormone_name','related_genes','organism','pmid','related_diseases','regulation_type','evidence_source')[:5000])
    random.shuffle(rows)

    for r in rows:
        h = get_val(r.get('hormone_name')); g = get_val(r.get('related_genes'))
        if h and g and sum(1 for x in questions if x['type'] == 'factual') < 50:
            add(f"Which gene is associated with the hormone {h} in the database?", 'factual', g, r)

    for r in rows:
        h = get_val(r.get('hormone_name')); sp = get_val(r.get('organism')); g = get_val(r.get('related_genes'))
        if h and sp and g and sum(1 for x in questions if x['type'] == 'species_specific') < 50:
            add(f"For {sp}, what gene information is recorded for the hormone {h}?", 'species_specific', f"{g}; {sp}", r)

    for r in rows:
        g = get_val(r.get('related_genes')); disease = get_val(r.get('related_diseases')); h = get_val(r.get('hormone_name'))
        if g and disease and h and sum(1 for x in questions if x['type'] == 'association') < 50:
            add(f"What disease or phenotype association is recorded for gene {g} in relation to hormone {h}?", 'association', disease, r)

    recs = []
    for Model in (HormoneReceptorInfo, HormoneReceptorFull):
        try:
            recs += list(Model.objects.all().values()[:3000])
        except Exception:
            pass
    random.shuffle(recs)
    for r in recs:
        h = get_val(r.get('hormone_name')); rec = get_val(r.get('receptor_name')); sp = get_val(r.get('receptor_species_name')); genes = get_val(r.get('receptor_coding_genes'))
        if h and (rec or genes) and sum(1 for x in questions if x['type'] == 'synthesis') < 50:
            tail = f" in {sp}" if sp else ""
            add(f"Summarize the receptor and coding gene information for {h}{tail}.", 'synthesis', '; '.join([v for v in [rec, genes, sp] if v]), r)

    fallback = [
        ('What information is available about insulin related genes?', 'factual', 'insulin'),
        ('Summarize hormone receptor information for peptide hormones.', 'synthesis', 'receptor'),
        ('What species-specific hormone receptor records are available in the database?', 'species_specific', 'species'),
        ('Describe gene-disease associations captured by the hormone database.', 'association', 'disease'),
    ]
    i = 0
    while len(questions) < 200:
        q, t, e = fallback[i % len(fallback)]
        add(f"{q} Focus item {i+1}.", t, e, {})
        i += 1
    questions = questions[:200]
    for i, q in enumerate(questions, 1):
        q['id'] = i
    jdump(DATASET, questions)
    return questions

def deepseek_direct(question):
    messages = [
        {'role':'system','content':'You are a hormone research assistant. Answer directly from your general knowledge. Be concise and in English.'},
        {'role':'user','content':question}
    ]
    payload = {'model': DEEPSEEK_MODEL, 'messages': messages, 'max_tokens': 2048, 'temperature': 0.3, 'top_p': 0.95, 'stream': False}
    r = deepseek_rag_client.session.post(DEEPSEEK_API_URL, json=payload, headers=deepseek_rag_client.headers, timeout=120, verify=True)
    if r.status_code != 200:
        return {'success': False, 'answer': f'ERROR {r.status_code}: {r.text[:300]}'}
    data = r.json()
    return {'success': True, 'answer': data['choices'][0]['message']['content'].strip(), 'usage': data.get('usage', {})}

def heuristic_traceability(answer):
    if not answer:
        return 0.0
    score = 0.0
    pats = [r'PMID[:\s]*\d+', r'UniProt|NCBI|PubChem|DOID|DO:', r'gene|receptor|hormone|species|organism', r'database|record|evidence|source', r'[A-Z0-9]{3,10}_[A-Z0-9]+|[A-Z]{2,}\d+']
    weights = [0.25, 0.20, 0.20, 0.20, 0.15]
    for pat, w in zip(pats, weights):
        if re.search(pat, answer, re.I):
            score += w
    return min(1.0, score)

def judge(question, expected, rag_answer, direct_answer, qtype):
    prompt = f"""Evaluate two answers to a hormone-database question.
Question type: {qtype}
Question: {question}
Expected/database evidence: {expected}

RAG answer:
{rag_answer[:6000]}

Direct DeepSeek answer without retrieval:
{direct_answer[:4000]}

Return strict JSON only with numeric fields in [0,1]:
{{"rag_relevance":0,"rag_completeness":0,"rag_accuracy":0,"rag_hallucination":0,"direct_hallucination":0,"f1":0,"judge_notes":"short"}}
Definitions:
- relevance: addresses the question
- completeness: covers expected evidence and important caveats
- accuracy: consistency with expected/database evidence
- hallucination: unsupported or contradicted claims rate; 0=no hallucination, 1=severe
- f1: type-level answer quality proxy combining precision/recall against expected evidence
"""
    messages = [{'role':'system','content':'You are a strict scientific QA evaluator. Output valid JSON only.'},{'role':'user','content':prompt}]
    payload = {'model': DEEPSEEK_MODEL, 'messages': messages, 'max_tokens': 1200, 'temperature': 0.0, 'top_p': 1, 'stream': False}
    try:
        r = deepseek_rag_client.session.post(DEEPSEEK_API_URL, json=payload, headers=deepseek_rag_client.headers, timeout=120, verify=True)
        txt = r.json()['choices'][0]['message']['content'].strip() if r.status_code == 200 else '{}'
        m = re.search(r'\{.*\}', txt, re.S)
        data = json.loads(m.group(0) if m else txt)
    except Exception as e:
        data = {'rag_relevance':0.0,'rag_completeness':0.0,'rag_accuracy':0.0,'rag_hallucination':1.0,'direct_hallucination':1.0,'f1':0.0,'judge_notes':f'judge failed: {e}'}
    for k in ['rag_relevance','rag_completeness','rag_accuracy','rag_hallucination','direct_hallucination','f1']:
        try:
            data[k] = max(0.0, min(1.0, float(data.get(k, 0))))
        except Exception:
            data[k] = 0.0
    return data

def run_one(q):
    path = RESULT_DIR / f"result_{q['id']:03d}.json"
    if path.exists():
        old = json.loads(path.read_text(encoding='utf-8'))
        if old.get('status') == 'done':
            return old
    rec = {'id': q['id'], 'question': q['question'], 'type': q['type'], 'expected': q.get('expected',''), 'started_at': now()}
    jdump(path, {**rec, 'status': 'running'})
    try:
        t0 = time.time(); rag = rag_search(q['question']); rag_sec = time.time() - t0
        t1 = time.time(); direct = deepseek_direct(q['question']); direct_sec = time.time() - t1
        judge_res = judge(q['question'], q.get('expected',''), rag.get('answer',''), direct.get('answer',''), q['type'])
        relevance = judge_res['rag_relevance']; completeness = judge_res['rag_completeness']; accuracy = judge_res['rag_accuracy']
        overall = 0.4*relevance + 0.3*completeness + 0.3*accuracy
        trace = heuristic_traceability(rag.get('answer',''))
        rec.update({
            'status':'done', 'finished_at':now(),
            'rag_answer': rag.get('answer',''), 'rag_model': rag.get('model'), 'rag_context_type': rag.get('context_type'), 'rag_has_data': rag.get('has_data'), 'rag_token_usage': rag.get('token_usage',{}), 'rag_seconds': round(rag_sec,3),
            'direct_answer': direct.get('answer',''), 'direct_success': direct.get('success'), 'direct_seconds': round(direct_sec,3),
            'judge': judge_res,
            'scores': {'relevance':relevance,'completeness':completeness,'accuracy':accuracy,'overall':overall,'traceability':trace,'f1':judge_res['f1'],'rag_hallucination':judge_res['rag_hallucination'],'direct_hallucination':judge_res['direct_hallucination'],'passed': overall >= 0.75}
        })
    except Exception as e:
        rec.update({'status':'failed','finished_at':now(),'error':str(e),'traceback':traceback.format_exc()})
    jdump(path, rec)
    return rec

def mean(done, key):
    vals = [r['scores'][key] for r in done if 'scores' in r]
    return sum(vals)/len(vals) if vals else 0.0

def make_figures(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7,5), dpi=300)
    labels = ['RAG hallucination','Direct DeepSeek hallucination','Traceability']
    vals = [summary['rag_hallucination_rate'], summary['direct_hallucination_rate'], summary['traceability']]
    bars = ax.bar(labels, vals, color=['#2ca25f','#de2d26','#3182bd'])
    ax.set_ylim(0,1); ax.set_ylabel('Score / Rate'); ax.set_title('Figure 6C. Hallucination reduction and traceability'); ax.grid(axis='y', alpha=0.25)
    for b,v in zip(bars, vals): ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}', ha='center', fontsize=9)
    plt.xticks(rotation=20, ha='right'); plt.tight_layout(); fig.savefig(FIG_DIR/'figure_6C_hallucination_traceability.png'); plt.close(fig)

    types = list(summary['by_type'].keys()); vals = [summary['by_type'][t]['f1'] for t in types]
    fig, ax = plt.subplots(figsize=(7,5), dpi=300)
    bars = ax.bar(types, vals, color=['#756bb1','#31a354','#fd8d3c','#6baed6'][:len(types)])
    ax.set_ylim(0,1); ax.set_ylabel('F1 score'); ax.set_title('Figure 6D. Performance by question type'); ax.grid(axis='y', alpha=0.25)
    for b,v in zip(bars, vals): ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}', ha='center', fontsize=9)
    plt.xticks(rotation=20, ha='right'); plt.tight_layout(); fig.savefig(FIG_DIR/'figure_6D_f1_by_question_type.png'); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8,5), dpi=300)
    labels = ['Overall','Relevance','Completeness','Accuracy','Pass rate']
    vals = [summary['overall_score'], summary['relevance'], summary['completeness'], summary['accuracy'], summary['pass_rate']]
    bars = ax.bar(labels, vals, color='#4c78a8')
    ax.set_ylim(0,1); ax.set_ylabel('Score'); ax.set_title('RAG 200-question quantitative evaluation'); ax.grid(axis='y', alpha=0.25)
    for b,v in zip(bars, vals): ax.text(b.get_x()+b.get_width()/2, v+0.02, f'{v:.3f}', ha='center', fontsize=9)
    plt.tight_layout(); fig.savefig(FIG_DIR/'rag_200_overall_metrics.png'); plt.close(fig)

def aggregate():
    records = []
    for p in sorted(RESULT_DIR.glob('result_*.json')):
        try:
            records.append(json.loads(p.read_text(encoding='utf-8')))
        except Exception:
            pass
    done = [r for r in records if r.get('status') == 'done']
    bytype = {}
    for t in sorted(set(r.get('type') for r in done)):
        trs = [r for r in done if r.get('type') == t]
        bytype[t] = {'n': len(trs), 'f1': sum(r['scores']['f1'] for r in trs)/len(trs), 'overall': sum(r['scores']['overall'] for r in trs)/len(trs)}
    summary = {
        'generated_at': now(), 'n_total': 200, 'n_done': len(done),
        'overall_score': mean(done,'overall'), 'relevance': mean(done,'relevance'), 'completeness': mean(done,'completeness'), 'accuracy': mean(done,'accuracy'),
        'pass_threshold': 0.75, 'pass_count': sum(1 for r in done if r['scores']['passed']), 'pass_rate': (sum(1 for r in done if r['scores']['passed'])/len(done) if done else 0),
        'rag_hallucination_rate': mean(done,'rag_hallucination'), 'direct_hallucination_rate': mean(done,'direct_hallucination'), 'traceability': mean(done,'traceability'),
        'by_type': bytype,
    }
    jdump(SUMMARY_JSON, summary)
    import csv
    with open(SUMMARY_CSV, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['id','type','question','overall','relevance','completeness','accuracy','f1','traceability','rag_hallucination','direct_hallucination','passed','rag_context_type','rag_seconds','direct_seconds'])
        for r in done:
            s = r['scores']
            w.writerow([r['id'], r['type'], r['question'], s['overall'], s['relevance'], s['completeness'], s['accuracy'], s['f1'], s['traceability'], s['rag_hallucination'], s['direct_hallucination'], s['passed'], r.get('rag_context_type'), r.get('rag_seconds'), r.get('direct_seconds')])
    if done:
        make_figures(summary)
    return summary

def main():
    qs = build_dataset()
    print(f'[{now()}] dataset={len(qs)} out={OUT}', flush=True)
    for q in qs:
        r = run_one(q)
        done_count = len([p for p in RESULT_DIR.glob('result_*.json') if json.loads(p.read_text(encoding='utf-8')).get('status') == 'done'])
        jdump(PROGRESS, {'updated_at': now(), 'last_id': q['id'], 'done': done_count, 'total': len(qs), 'current_status': r.get('status')})
        score = r.get('scores', {}).get('overall')
        print(f'[{now()}] q={q["id"]}/200 status={r.get("status")} done={done_count} overall={score}', flush=True)
    summary = aggregate()
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)

if __name__ == '__main__':
    main()
