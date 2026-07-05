# raredx — Análisis de VCF con IA para diagnóstico de enfermedades raras

Pipeline moderno que anota un VCF y prioriza variantes candidatas para diagnóstico de
enfermedades raras, combinando **cuatro capas de evidencia + dos capas de IA**, con
priorización dirigida por fenotipo al estilo de Exomiser.

## Arquitectura

```
  Nota clínica (texto libre)
        |  [IA-1] LLM extrae fenotipos -> OLS4 los normaliza a HPO
        v
  Perfil HPO del paciente --------------+
                                        |
  VCF --> anotacion por variante        |  MOTOR DE FENOTIPO
        |- Ensembl VEP (consecuencia)   |  gen->enfermedad->HPO (Open Targets)
        |- gnomAD r4 (frecuencia, pLI)  |  + expansion ontologica
        |- ClinVar (significancia, *)   |  -> pheno_score
        |- [IA-2] ESM-2 (missense LLR)  |
        +-> MOTOR DE VARIANTE (ACMG)    |
              -> variant_score ---------+
                                        v
                     combined = 0.55*variante + 0.45*fenotipo
                                        v
                        Ranking + informe HTML + CSV
```

## Las dos capas de IA (novedad)

### IA-1 - Extraccion de fenotipo HPO desde nota clinica (LLM)
En vez de introducir codigos HPO a mano, **pegas la nota clinica en lenguaje natural** y un
LLM (Claude) extrae los fenotipos presentes, que OLS4 aterriza a terminos HPO oficiales.
Captura signos sutiles: *"el nino sabe salado"* -> `HP:0012236` (cloruro en sudor elevado,
patognomonico de fibrosis quistica), *"acropaquias"* -> `HP:0100759` (dedos en palillo de tambor).

```bash
python raredx_pipeline.py input.vcf --clinical-note examples/clinical_note_es.txt --out-prefix out/p
```

### IA-2 - Patogenicidad de missense con ESM-2 (modelo de lenguaje de proteinas)
Reemplaza/complementa SIFT/PolyPhen (predictores de ~2010) con **ESM-2** (Meta AI), que
puntua cada cambio de aminoacido por su verosimilitud evolutiva (masked-marginal
log-likelihood ratio). Un LLR muy negativo => mutacion improbable => probablemente deleterea.
Alimenta las reglas **PP3/BP4** del motor ACMG. Es el principio detras de AlphaMissense.

```bash
python raredx_pipeline.py input.vcf --esm --hpo "HP:0002205,HP:0001738" --out-prefix out/p
```

## Motor de variante (ACMG-lite)

| Capa | Fuente | Aporta |
|------|--------|--------|
| Consecuencia funcional | Ensembl VEP | missense/stop/frameshift + SIFT/PolyPhen |
| Frecuencia poblacional | gnomAD r4 | AF (BA1/BS1: comun => benigna) |
| Restriccion genica | gnomAD | pLI, LOEUF (PVS1 en genes intolerantes a LoF) |
| Significancia clinica | ClinVar | veredicto + estrellas de oro |
| **IA missense** | **ESM-2** | **LLR evolutivo -> PP3/BP4** |

## Efecto demostrado

Con un perfil de fibrosis quistica extraido de la nota clinica (10 terminos HPO), **CFTR
sube del puesto #3-4 (solo variante) al #1**; BRCA1 -mayor score de variante pero sin
relacion fenotipica- baja al #2. Ademas, ESM-2 baja el score de BRCA1 A566V (LLR +0.21,
tolerada), refinando un VUS. El diseno multicapa evita falsos positivos: APOE R176C sale
deleterea por ESM, pero la frecuencia (BA1) y el fenotipo la mantienen correctamente benigna.

## Uso completo

```bash
python raredx_pipeline.py input.vcf \
       --clinical-note nota.txt \      # IA-1: HPO desde texto (o --hpo "HP:...")
       --esm \                          # IA-2: ESM-2 en missense
       --out-prefix salida/paciente \
       --email tu@institucion.org
```

**Dependencias:** `requests` (base); `torch fair-esm` (para `--esm`); `anthropic` +
`ANTHROPIC_API_KEY` (para `--clinical-note`). Sin GPU, ESM-2 usa el modelo 8M en CPU.

## Salidas
- `<prefix>_report.html` - informe clinico: nota->HPO, ranking con LLR de ESM-2, fichas de evidencia
- `<prefix>_annotated.csv` - tabla completa con scores de variante, ESM-2, fenotipo y combinado

## Nota importante
Sistema de **apoyo a la decision, no un diagnostico**. Todo hallazgo debe confirmarse por
metodo ortogonal e interpretarse por un genetista clinico en el contexto completo del paciente.
Los scores de ESM-2 con el modelo 8M son una prueba de concepto; para uso serio, emplear
ESM-2 650M/3B (GPU) o AlphaMissense.
