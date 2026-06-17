# PROJECT_MEMORY — Animal Hormone DB / RAG Assistant

## Project identity
- Server: `43.99.62.219`
- Production path: `/www/wwwroot/default/animal_hormone`
- Main service: Django + Gunicorn, public Animal Hormone DB site.
- RAG implementation files:
  - `hormone_app/utils/rag_deepseek.py`
  - `hormone_app/utils/rag_retriever.py`
- DeepSeek model configured in code: `deepseek-v4-pro`.

## 2026-05-31 — RAG 200-question quantitative evaluation
ouev5 requested a quantitative evaluation matching the manuscript paragraph:
- 200-question test set
- save one result after each question
- compare RAG vs direct DeepSeek generation
- compute relevance, completeness, accuracy, overall score, pass rate, hallucination rates, traceability and F1 by question type
- generate figures equivalent to Figure 6C and Figure 6D

Implemented evaluation workspace:
- Directory: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/`
- Script: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/run_rag_200_eval.py`
- Dataset path: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/rag_200_questions.json`
- Per-question results: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/results/result_XXX.json`
- Progress file: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/progress.json`
- Logs: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/logs/`
- Summary outputs after completion:
  - `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/rag_200_summary.json`
  - `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/rag_200_results.csv`
- Figure outputs after completion:
  - `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/figures/figure_6C_hallucination_traceability.png`
  - `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/figures/figure_6D_f1_by_question_type.png`
  - `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval/figures/rag_200_overall_metrics.png`

Current verified status during setup:
- Evaluation started successfully.
- Vector index loaded from precomputed embeddings: 18,102 records, 384 dimensions.
- At last check: at least `2/200` questions completed; question 3 had started.
- Example progress observed: `q=2/200 status=done done=2 overall=1.0`.
- Important: this is a long-running evaluation; do not claim final manuscript numbers until all 200 result JSON files and final summary/figures exist.

How to check status:
```bash
cd /www/wwwroot/default/animal_hormone
pgrep -af "run_rag_200_eval"
cat evaluation/rag_200_eval/progress.json
find evaluation/rag_200_eval/results -name "result_*.json" | wc -l
ls -lt evaluation/rag_200_eval/logs | head
ls -l evaluation/rag_200_eval/rag_200_summary.json evaluation/rag_200_eval/figures/*.png 2>/dev/null
```

How to resume safely:
```bash
cd /www/wwwroot/default/animal_hormone
nohup venv/bin/python evaluation/rag_200_eval/run_rag_200_eval.py \
  > evaluation/rag_200_eval/logs/run_$(date +%Y%m%d_%H%M%S).log 2>&1 &
```
The script is resumable: existing `status=done` result files are skipped.

Caution:
- The first implementation accidentally ran once via `--help` because the script did not parse CLI args; this still produced valid per-question files, but pid writing was unreliable. Use status files/logs/results count as source of truth.
- The evaluation uses DeepSeek for RAG answer, direct baseline, and judge scoring. Report actual computed results only after completion.

## Update — first 200-question run completed, but metrics are not manuscript-ready
- Completion verified: `progress.json` shows `done=200/200` at `2026-05-31T08:08:22`.
- Summary and figures were generated under `evaluation/rag_200_eval/`.
- Actual first-run metrics were very low: overall `0.072`, pass rate `7/200`, RAG hallucination `0.96`, direct hallucination `0.985`.
- These numbers likely reflect a flawed/misaligned evaluation design rather than final system quality, because expected answers were automatically derived from sparse fields and the judge penalized heavily.
- Do not cite this run as the paper's final quantitative evaluation. Revise evaluation protocol and rerun before writing manuscript results.

## RAG 200-question evaluation v3 started — field-anchored protocol
- ouev5 approved fixing the evaluation protocol rather than fabricating high metrics.
- Inspected low-scoring v1 cases and found systematic false negatives: RAG answers often contained the expected field (e.g. AR), but the LLM judge scored 0; original questions were too broad (`which gene is associated`) and gold answers used only one row/field despite multiple valid database associations.
- Created v2 deterministic recomputation from saved answers under `evaluation/rag_200_eval_v2/`; it improved overall from `0.072` to `0.831` but still exposed weak association/gold mismatch.
- Created v3 runner: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval_v3_runner.py`.
- v3 output directory: `/www/wwwroot/default/animal_hormone/evaluation/rag_200_eval_v3/`.
- v3 design:
  - regenerates 200 explicit field-anchored questions;
  - questions specify exact database fields such as `related_genes`, `related disease/phenotype`, `receptor_name`, `receptor_coding_genes`;
  - runs RAG and direct DeepSeek baseline;
  - saves one JSON after each question under `evaluation/rag_200_eval_v3/results/result_XXX.json`;
  - continuously writes `progress.json`, `rag_200_summary_v3.json`, `rag_200_results_v3.csv`, and figures.
- v3 started in background with PID recorded in `evaluation/rag_200_eval_v3/runner.pid`.
- Early sanity check after 2/200:
  - overall_score: `0.93346`
  - relevance: `0.89637`
  - completeness: `0.93204`
  - accuracy: `0.98433`
  - pass_count/pass_rate: `2/2 = 1.0`
  - rag_hallucination_rate: `0.015`
  - direct_hallucination_rate: `0.88167`
  - traceability: `0.85`
- Do not report final v3 metrics until `progress.json` shows `done=200` and final summary/figures are generated.

## RAG 200-question evaluation v3 completed
- Final v3 progress checked on 2026-05-31 12:49 CST.
- `evaluation/rag_200_eval_v3/progress.json` shows `done=200/200`, `last_id=200`, updated at `2026-05-31T12:30:59`.
- Result count: 200 `result_*.json` files.
- Runner process ended after completion.
- Final v3 summary:
  - overall_score: `0.7349569542124543`
  - relevance: `0.9268623511904762`
  - completeness: `0.6420021153846155`
  - accuracy: `0.5720379304029304`
  - pass_count/pass_rate: `100/200 = 0.5`
  - rag_hallucination_rate: `0.4019559065934066`
  - direct_hallucination_rate: `0.7369571503496504`
  - traceability: `0.74725`
- By type:
  - factual: n=50, F1=`0.7606061172620572`, overall=`0.841334336996337`
  - species_specific: n=50, F1=`0.6069298795868798`, overall=`0.7489086855366855`
  - association: n=50, F1=`0.20077696540098244`, overall=`0.5287347383727384`
  - synthesis: n=50, F1=`0.7797449095291421`, overall=`0.820850055944056`
- Interpretation: v3 completed and improved evaluation design, but final metrics still do not support ouev5's desired high-performance manuscript paragraph. Main weak area remains association questions. Next step should be inspecting association failures and improving retrieval/field matching or revising association gold protocol before manuscript use.
