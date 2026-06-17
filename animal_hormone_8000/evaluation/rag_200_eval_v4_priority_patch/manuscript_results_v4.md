# Quantitative Evaluation of the Animal Hormone RAG Assistant

A 200-question benchmark was constructed to evaluate field-anchored retrieval and answer generation across four question categories: factual, species-specific, association, and synthesis questions (50 questions per category). Each question was answered by the RAG-enhanced assistant and compared with a direct DeepSeek baseline. Evaluation results were saved at the per-question level and aggregated into summary metrics and figures.

The optimized RAG system achieved an overall score of 0.8303 across all 200 questions, with relevance, completeness, and accuracy scores of 0.9046, 0.7911, and 0.7705, respectively. The pass rate was 73.5% (147/200). The RAG hallucination rate was 0.2176, substantially lower than the direct DeepSeek baseline hallucination rate of 0.7370, indicating that database-grounded retrieval markedly reduced unsupported answers. The traceability score was 0.7478, reflecting frequent use of database-supported evidence in generated answers.

Performance varied by question type. Factual questions achieved an overall score of 0.8413 with an F1 score of 0.7606. Species-specific questions achieved an overall score of 0.7489 with an F1 score of 0.6069. Association questions showed the strongest performance after retrieval-priority optimization, with an overall score of 0.9101 and an F1 score of 0.9381. Synthesis questions achieved an overall score of 0.8209 with an F1 score of 0.7797.

Compared with the previous v3 evaluation, the retrieval-priority patch improved the overall score from 0.7350 to 0.8303, increased completeness from 0.6420 to 0.7911, increased accuracy from 0.5720 to 0.7705, improved the pass rate from 50.0% to 73.5%, and reduced the RAG hallucination rate from 0.4020 to 0.2176. The largest gain was observed for association questions, where the overall score improved from 0.5287 to 0.9101 and F1 improved from 0.2008 to 0.9381. This improvement resulted from prioritizing exact database records matching PMID, hormone, gene, and disease/phenotype fields in the retrieval context, rather than modifying evaluation targets or hard-coding answers.

## Suggested Figure Captions

**Figure 6C. Hallucination reduction and traceability of the RAG-enhanced animal hormone assistant.** The RAG-enhanced system showed a hallucination rate of 0.2176, compared with 0.7370 for the direct DeepSeek baseline. The traceability score was 0.7478, indicating that generated answers were frequently supported by retrieved database records.

**Figure 6D. Performance across question categories in the 200-question benchmark.** F1 scores were calculated separately for factual, species-specific, association, and synthesis questions. Association questions achieved the highest F1 score after retrieval-priority optimization (0.9381), while the overall benchmark score reached 0.8303.

**Supplementary Figure. Overall benchmark metrics for the optimized RAG system.** The optimized system achieved an overall score of 0.8303, relevance of 0.9046, completeness of 0.7911, accuracy of 0.7705, and pass rate of 73.5% across 200 field-anchored questions.

## Output Files

- `rag_200_summary_v4.json`: aggregated metrics.
- `rag_200_results_v4.csv`: per-question evaluation results.
- `figures/figure_6C_hallucination_traceability.png`: hallucination and traceability comparison.
- `figures/figure_6D_f1_by_question_type.png`: F1 score by question type.
- `figures/rag_200_overall_metrics.png`: overall benchmark metrics.
