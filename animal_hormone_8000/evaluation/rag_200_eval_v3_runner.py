#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""RAG 200-question v3 evaluation with field-anchored gold answers.
- Generates explicit questions tied to specific DB fields.
- Saves one result after each question.
- Runs RAG answer and direct DeepSeek baseline.
- Uses deterministic gold-field scoring and continuously updates summary/figures.
"""
import os, sys, json, csv, time, random, re, math, traceback
from pathlib import Path
from datetime import datetime

BASE = Path('/www/wwwroot/default/animal_hormone')
OUT = BASE / 'evaluation' / 'rag_200_eval_v3'
RESULT_DIR = OUT / 'results'
FIG = OUT / 'figures'
LOG = OUT / 'logs'
DATASET = OUT / 'rag_200_questions_v3.json'
SUMMARY = OUT / 'rag_200_summary_v3.json'
CSV = OUT / 'rag_200_results_v3.csv'
PROGRESS = OUT / 'progress.json'
for d in [OUT, RESULT_DIR, FIG, LOG]: d.mkdir(parents=True, exist_ok=True)

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'animal_hormone.settings')
sys.path.insert(0, str(BASE))
import django
django.setup()
from hormone_app.models import HormoneRelatedGene, HormoneReceptorInfo, HormoneReceptorFull
from hormone_app.utils.rag_deepseek import rag_search, deepseek_rag_client, DEEPSEEK_API_URL, DEEPSEEK_MODEL

random.seed(20260531)

def now(): return datetime.now().isoformat(timespec='seconds')
def jdump(path, obj):
    tmp=Path(str(path)+'.tmp'); tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2), encoding='utf-8'); tmp.replace(path)
def val(x): return str(x or '').strip()
def norm(s): return re.sub(r'[^a-z0-9]+',' ',str(s or '').lower()).strip()
def words(s): return [w for w in norm(s).split() if len(w)>=2]
def clamp(x): return max(0.0, min(1.0, float(x)))
def has_term(ans, term):
    ans=ans or ''; term=str(term or '').strip()
    if not term: return False
    if re.fullmatch(r'[A-Za-z0-9_.-]{2,30}', term):
        return re.search(r'\b'+re.escape(term)+r'\b', ans, re.I) is not None
    a=norm(ans); t=norm(term)
    if re.search(r'(^| )'+re.escape(t)+r'( |$)', a): return True
    tw=words(term)
    if not tw: return False
    present=sum(1 for w in tw if re.search(r'(^| )'+re.escape(w)+r'( |$)', a))
    return present >= max(1, math.ceil(len(tw)*0.70))
def split_terms(s):
    terms=[]
    for p in re.split(r'[;,|]+|\band\b', str(s or ''), flags=re.I):
        p=p.strip()
        if p and len(p)<=100: terms.append(p)
    out=[]; seen=set()
    for t in terms:
        k=norm(t)
        if k and k not in seen:
            out.append(t); seen.add(k)
    return out[:10]
def coverage(ans, expected_terms):
    if not expected_terms: return 0.5
    return sum(1 for t in expected_terms if has_term(ans,t))/len(expected_terms)
def focus_cov(q, ans):
    stop=set('according animal hormone database what which is are the field value values listed recorded for in of and or to with record pmid answer provide based on'.split())
    qs=[w for w in words(q) if w not in stop and len(w)>2]
    if not qs: return 0.8 if ans else 0
    a=norm(ans)
    return clamp(sum(1 for w in qs[:15] if re.search(r'(^| )'+re.escape(w)+r'( |$)', a))/min(len(qs),15))
def trace(ans):
    pats=[r'PMID[:\s]*\d+', r'UniProt|NCBI|PubChem|DOID|DO:', r'database|record|field|listed|provided context', r'gene|receptor|hormone|species|organism', r'\b[A-Z][A-Za-z0-9_.-]{2,12}\b']
    weights=[.25,.15,.25,.2,.15]
    return clamp(sum(w for pat,w in zip(pats,weights) if re.search(pat, ans or '', re.I)))
def quality(ans):
    if not ans: return 0
    low=ans.lower()
    if low.startswith('error') or 'request timeout' in low[:120]: return 0
    return clamp(0.75*clamp(len(ans)/260) + (0.15 if any(c in ans for c in [':','-','\n','•']) else 0) + (0.10 if 'database' in low or 'record' in low else 0))

def build_dataset():
    if DATASET.exists(): return json.loads(DATASET.read_text(encoding='utf-8'))
    qs=[]; seen=set()
    def add(qtype, question, expected, evidence):
        if not expected: return
        key=norm(question)
        if key in seen: return
        seen.add(key)
        qs.append({'id':len(qs)+1,'type':qtype,'question':question,'expected':expected,'expected_terms':split_terms(expected),'evidence':evidence})
    rows=list(HormoneRelatedGene.objects.exclude(hormone_name__isnull=True).exclude(related_genes__isnull=True).values('id','hormone_name','related_genes','organism','pmid','related_diseases','regulation_type','evidence_source')[:10000])
    random.shuffle(rows)
    for r in rows:
        if sum(1 for x in qs if x['type']=='factual')>=50: break
        h=val(r.get('hormone_name')); g=val(r.get('related_genes')); org=val(r.get('organism')); pmid=val(r.get('pmid'))
        if h and g:
            tail=f" in {org}" if org else ''
            pm=f" with PMID {pmid}" if pmid else ''
            add('factual', f"According to the Animal Hormone Database, what is the related_genes field value for hormone {h}{tail}{pm}?", g, r)
    for r in rows:
        if sum(1 for x in qs if x['type']=='species_specific')>=50: break
        h=val(r.get('hormone_name')); g=val(r.get('related_genes')); org=val(r.get('organism'))
        if h and g and org:
            add('species_specific', f"For organism {org}, which related gene is recorded for hormone {h} in the Animal Hormone Database related_genes field?", g, r)
    for r in rows:
        if sum(1 for x in qs if x['type']=='association')>=50: break
        h=val(r.get('hormone_name')); g=val(r.get('related_genes')); dis=val(r.get('related_diseases')); pmid=val(r.get('pmid'))
        if h and g and dis and dis.lower() not in {'none','nan','null','not available',''}:
            pm=f" with PMID {pmid}" if pmid else ''
            add('association', f"According to the database record, what related disease or phenotype is listed for hormone {h} and gene {g}{pm}?", dis, r)
    recs=[]
    for Model in (HormoneReceptorInfo, HormoneReceptorFull):
        try: recs += list(Model.objects.all().values()[:6000])
        except Exception: pass
    random.shuffle(recs)
    for r in recs:
        if sum(1 for x in qs if x['type']=='synthesis')>=50: break
        h=val(r.get('hormone_name')); org=val(r.get('receptor_species_name') or r.get('hormone_species_name')); rec=val(r.get('receptor_name')); gene=val(r.get('receptor_coding_genes'))
        expected='; '.join([x for x in [rec,gene] if x])
        if h and expected:
            tail=f" in {org}" if org else ''
            add('synthesis', f"Summarize the receptor_name and receptor_coding_genes fields for hormone {h}{tail} according to the Animal Hormone Database.", expected, r)
    qs=qs[:200]
    for i,q in enumerate(qs,1): q['id']=i
    if len(qs)!=200: raise RuntimeError(f'Only built {len(qs)} questions')
    jdump(DATASET, qs)
    return qs

def direct_answer(q):
    messages=[{'role':'system','content':'You are a hormone research assistant. Answer directly from general knowledge only, without database retrieval. English only.'},{'role':'user','content':q}]
    payload={'model':DEEPSEEK_MODEL,'messages':messages,'max_tokens':1536,'temperature':0.3,'top_p':0.95,'stream':False}
    r=deepseek_rag_client.session.post(DEEPSEEK_API_URL,json=payload,headers=deepseek_rag_client.headers,timeout=120,verify=True)
    if r.status_code!=200: return {'success':False,'answer':f'ERROR {r.status_code}: {r.text[:250]}'}
    data=r.json(); return {'success':True,'answer':data['choices'][0]['message']['content'].strip(),'usage':data.get('usage',{})}

def score(q, rag_ans, direct_ans):
    terms=q.get('expected_terms') or split_terms(q.get('expected',''))
    cov=coverage(rag_ans, terms); dcov=coverage(direct_ans, terms); foc=focus_cov(q['question'], rag_ans); tr=trace(rag_ans); qual=quality(rag_ans)
    relevance=clamp(.58*foc+.27*qual+.15)
    accuracy=clamp(.86*cov+.08*foc+.06*tr)
    completeness=clamp(.62*cov+.20*qual+.18*tr)
    overall=clamp(.4*relevance+.3*completeness+.3*accuracy)
    f1=clamp(2*accuracy*completeness/(accuracy+completeness)) if accuracy+completeness else 0
    rag_h=clamp(.82*(1-cov)+.18*(1-foc))
    direct_h=clamp(.86*(1-dcov)+.14*(1-focus_cov(q['question'], direct_ans)))
    return {'expected_coverage':cov,'direct_expected_coverage':dcov,'relevance':relevance,'completeness':completeness,'accuracy':accuracy,'overall':overall,'traceability':tr,'f1':f1,'rag_hallucination':rag_h,'direct_hallucination':direct_h,'passed':overall>=.75,'focus_coverage':foc,'answer_quality':qual}

def make_figures(summary):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,5),dpi=300)
    labels=['RAG hallucination','Direct DeepSeek hallucination','Traceability']; vals=[summary['rag_hallucination_rate'],summary['direct_hallucination_rate'],summary['traceability']]
    bars=ax.bar(labels,vals,color=['#2ca25f','#de2d26','#3182bd']); ax.set_ylim(0,1); ax.set_ylabel('Score / Rate'); ax.set_title('Figure 6C. Hallucination reduction and traceability'); ax.grid(axis='y',alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+.02,f'{v:.3f}',ha='center',fontsize=9)
    plt.xticks(rotation=18,ha='right'); plt.tight_layout(); fig.savefig(FIG/'figure_6C_hallucination_traceability.png'); plt.close(fig)
    types=list(summary['by_type']); vals=[summary['by_type'][t]['f1'] for t in types]
    fig,ax=plt.subplots(figsize=(7,5),dpi=300); bars=ax.bar(types,vals,color=['#756bb1','#31a354','#fd8d3c','#6baed6'][:len(types)]); ax.set_ylim(0,1); ax.set_ylabel('F1 score'); ax.set_title('Figure 6D. Performance by question type'); ax.grid(axis='y',alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+.02,f'{v:.3f}',ha='center',fontsize=9)
    plt.xticks(rotation=20,ha='right'); plt.tight_layout(); fig.savefig(FIG/'figure_6D_f1_by_question_type.png'); plt.close(fig)
    labels=['Overall','Relevance','Completeness','Accuracy','Pass rate']; vals=[summary['overall_score'],summary['relevance'],summary['completeness'],summary['accuracy'],summary['pass_rate']]
    fig,ax=plt.subplots(figsize=(8,5),dpi=300); bars=ax.bar(labels,vals,color='#4c78a8'); ax.set_ylim(0,1); ax.set_ylabel('Score'); ax.set_title('RAG 200-question field-anchored evaluation'); ax.grid(axis='y',alpha=.25)
    for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v+.02,f'{v:.3f}',ha='center',fontsize=9)
    plt.tight_layout(); fig.savefig(FIG/'rag_200_overall_metrics.png'); plt.close(fig)

def aggregate():
    recs=[]
    for p in sorted(RESULT_DIR.glob('result_*.json')):
        try:
            d=json.loads(p.read_text(encoding='utf-8'))
            if d.get('status')=='done': recs.append(d)
        except Exception: pass
    if not recs: return None
    def mean(k): return sum(r['scores'][k] for r in recs)/len(recs)
    by={}
    for t in ['factual','species_specific','association','synthesis']:
        xs=[r for r in recs if r['type']==t]
        if xs: by[t]={'n':len(xs),'f1':sum(r['scores']['f1'] for r in xs)/len(xs),'overall':sum(r['scores']['overall'] for r in xs)/len(xs)}
    def pval(a,n,b,m):
        if not n or not m: return None
        p=(a+b)/(n+m); se=math.sqrt(p*(1-p)*(1/n+1/m)) if 0<p<1 else 0
        return 1.0 if se==0 else math.erfc(abs((a/n-b/m)/se)/math.sqrt(2))
    rf=sum(1 for r in recs if r['scores']['rag_hallucination']>=.5); df=sum(1 for r in recs if r['scores']['direct_hallucination']>=.5)
    summary={'generated_at':now(),'rubric':'v3 field-anchored gold answers; one JSON per question','n_total':200,'n_done':len(recs),'overall_score':mean('overall'),'relevance':mean('relevance'),'completeness':mean('completeness'),'accuracy':mean('accuracy'),'pass_threshold':.75,'pass_count':sum(1 for r in recs if r['scores']['passed']),'pass_rate':sum(1 for r in recs if r['scores']['passed'])/len(recs),'rag_hallucination_rate':mean('rag_hallucination'),'direct_hallucination_rate':mean('direct_hallucination'),'rag_hallucination_flag_rate':rf/len(recs),'direct_hallucination_flag_rate':df/len(recs),'hallucination_flag_p_value':pval(rf,len(recs),df,len(recs)),'traceability':mean('traceability'),'by_type':by}
    jdump(SUMMARY, summary)
    with open(CSV,'w',newline='',encoding='utf-8') as f:
        w=csv.writer(f); w.writerow(['id','type','question','expected','overall','relevance','completeness','accuracy','f1','traceability','rag_hallucination','direct_hallucination','passed','expected_coverage','direct_expected_coverage','rag_seconds','direct_seconds'])
        for r in recs:
            s=r['scores']; w.writerow([r['id'],r['type'],r['question'],r['expected'],s['overall'],s['relevance'],s['completeness'],s['accuracy'],s['f1'],s['traceability'],s['rag_hallucination'],s['direct_hallucination'],s['passed'],s['expected_coverage'],s['direct_expected_coverage'],r.get('rag_seconds'),r.get('direct_seconds')])
    make_figures(summary)
    return summary

def run_one(q):
    p=RESULT_DIR/f"result_{q['id']:03d}.json"
    if p.exists():
        d=json.loads(p.read_text(encoding='utf-8'))
        if d.get('status')=='done': return d
    rec={**q,'status':'running','started_at':now()}; jdump(p,rec)
    try:
        t0=time.time(); rag=rag_search(q['question']); rs=time.time()-t0
        t1=time.time(); direct=direct_answer(q['question']); ds=time.time()-t1
        rag_ans=rag.get('answer',''); direct_ans=direct.get('answer','')
        rec.update({'status':'done','finished_at':now(),'rag_answer':rag_ans,'direct_answer':direct_ans,'rag_context_type':rag.get('context_type'),'rag_has_data':rag.get('has_data'),'rag_token_usage':rag.get('token_usage',{}),'direct_success':direct.get('success'),'direct_usage':direct.get('usage',{}),'rag_seconds':round(rs,3),'direct_seconds':round(ds,3),'scores':score(q,rag_ans,direct_ans)})
    except Exception as e:
        rec.update({'status':'failed','finished_at':now(),'error':str(e),'traceback':traceback.format_exc()})
    jdump(p,rec); return rec

def main():
    qs=build_dataset(); print(f'[{now()}] v3 dataset={len(qs)} out={OUT}', flush=True)
    for q in qs:
        r=run_one(q); done=len([p for p in RESULT_DIR.glob('result_*.json') if json.loads(p.read_text(encoding='utf-8')).get('status')=='done'])
        summary=aggregate()
        jdump(PROGRESS, {'updated_at':now(),'last_id':q['id'],'done':done,'total':len(qs),'current_status':r.get('status'),'partial_summary':summary})
        print(f'[{now()}] q={q["id"]}/200 status={r.get("status")} done={done} overall={r.get("scores",{}).get("overall")}', flush=True)
    print(json.dumps(aggregate(), ensure_ascii=False, indent=2), flush=True)

if __name__=='__main__': main()
