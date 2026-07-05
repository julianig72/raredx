# raredx — Sistema de análisis de VCF con priorización por fenotipo (HPO)

Pipeline moderno y transparente que anota un VCF y prioriza variantes candidatas
para diagnóstico de enfermedades raras, combinando evidencia **a nivel de variante**
(estilo ACMG/AMP) con **coincidencia fenotípica HPO** paciente↔gen, al estilo de Exomiser.

## Arquitectura de dos motores

### 1. Motor de variante (ACMG-lite)
Para cada variante:

| Capa | Fuente | Campos |
|------|--------|--------|
| Consecuencia funcional | Ensembl VEP | consecuencia más severa, impacto, SIFT/PolyPhen |
| Frecuencia poblacional | gnomAD r4 | AF exoma/genoma, homocigotos |
| Restricción génica | gnomAD constraint | pLI, LOEUF |
| Significancia clínica | ClinVar (NCBI) | clasificación germinal, estrellas de oro, condición |

Reglas: **BA1/BS1** (AF ≥5%/≥1% → benigno), **PM2** (ausente/raro), **PVS1** (variante nula
en gen intolerante a LoF), **PP3/BP4** (in silico), y evidencia ClinVar ponderada por estrellas.

### 2. Motor de fenotipo (HPO)
- Los términos HPO del paciente (IDs `HP:xxxxxxx` o texto libre) se **resuelven vía OLS4**
  y se **expanden por la ontología** (ancestros) para un emparejamiento generoso.
- Para cada gen candidato, se recuperan las **enfermedades asociadas (Open Targets)** y sus
  **fenotipos HPO**, y se cruzan con el perfil del paciente.
- `pheno_score` = (coincidencias directas ×1.0 + coincidencias por ancestro ×0.4) / nº términos del paciente.

### Ranking combinado
`combined = 0.55·score_variante + 0.45·score_fenotipo`  (penalización ×0.5 si falla QC).

**Efecto demostrado:** con un perfil fenotípico tipo fibrosis quística (7 términos HPO:
infecciones respiratorias recurrentes, bronquiectasias, insuficiencia pancreática exocrina,
íleo meconial, cloruro en sudor elevado, retraso del crecimiento, tos crónica), **CFTR sube
del puesto #3–4 (solo variante, empatado con MYOC) al #1**, mientras que BRCA1 —con mayor score de variante pero
sin relación fenotípica con este paciente— baja al #2. Esta es la ganancia de rendimiento
diagnóstico que aporta el fenotipo.

## Uso del script

```bash
# Solo anotación + ranking ACMG:
python raredx_pipeline.py input.vcf --sample PATIENT_001 --out-prefix out/patient

# Con priorización por fenotipo (IDs HPO o texto libre, coma-separados o @archivo):
python raredx_pipeline.py input.vcf \
       --hpo "HP:0002205,HP:0001738,Bronchiectasis" \
       --out-prefix out/patient --email tu@institucion.org
```

El script es **autónomo y funcional de extremo a extremo**: anota contra las APIs públicas
en vivo (Ensembl VEP · gnomAD GraphQL · ClinVar E-utilities · Open Targets · HPO/OLS4),
sin necesidad de clave de API. Solo requiere `requests`.

## Salidas

- `patient_001_report.html` — informe clínico con perfil HPO, tabla de ranking combinado y fichas de evidencia por variante
- `patient_001_annotated.csv` — tabla completa: anotación + scores de variante, fenotipo y combinado
- `patient_001.vcf` — VCF de demostración (10 variantes reales GRCh38)
- `raredx_pipeline.py` — pipeline reutilizable (corre en cualquier VCF)

## Nota importante

Sistema de **apoyo a la decisión**, no un diagnóstico. Todo hallazgo reportable debe
confirmarse por un método ortogonal e interpretarse por un genetista clínico en el contexto
fenotípico completo del paciente. Los términos HPO deben capturarse cuidadosamente en la
consulta — la calidad del perfil fenotípico determina la calidad de la priorización.
