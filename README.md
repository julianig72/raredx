# raredx — Análisis de VCF con IA para diagnóstico de enfermedades raras

Pipeline moderno que anota un VCF y prioriza variantes candidatas para diagnóstico de
enfermedades raras, combinando **cuatro capas de evidencia + tres capas de IA**, con
priorización dirigida por fenotipo al estilo de Exomiser. Soporta **GRCh38 y GRCh37/hg19**
(este último con liftover automático a GRCh38 donde hace falta).

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
        |- [IA-3] AlphaMissense (0-1)   |
        +-> MOTOR DE VARIANTE (ACMG)    |
              -> variant_score ---------+
                                        v
                     combined = 0.55*variante + 0.45*fenotipo
                                        v
                        Ranking + informe HTML + CSV
```

## Las tres capas de IA (novedad)

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

### IA-3 - Patogenicidad de missense con AlphaMissense (DeepMind)
**AlphaMissense** (Cheng et al. 2023, *Science*) es un predictor calibrado clínicamente que
clasifica cada missense humana como `likely_pathogenic`/`ambiguous`/`likely_benign` con una
puntuación 0-1. DeepMind **no liberó el modelo ejecutable** (licencia no comercial), solo las
puntuaciones precalculadas de ~71 M variantes; el pipeline las consulta vía **Ensembl VEP**
(flag `AlphaMissense=1`). Como AlphaMissense es **solo GRCh38**, un VCF GRCh37 se **lifta a
GRCh38** por variante antes de consultar. Alimenta las reglas **PP3/BP4** (etiquetas `PP3_AM`/`BP4_AM`)
y pesa algo más que ESM-2 8M en el score de variante, por ser un predictor clínico especializado.
No requiere GPU ni descargar pesos.

```bash
python raredx_pipeline.py input.vcf --alphamissense --assembly GRCh37 --hpo "HP:0001250" --out-prefix out/p
```

## Motor de variante (ACMG-lite)

| Capa | Fuente | Aporta |
|------|--------|--------|
| Consecuencia funcional | Ensembl VEP | missense/stop/frameshift + SIFT/PolyPhen |
| Frecuencia poblacional | gnomAD r4 | AF (BA1/BS1: comun => benigna) |
| Restriccion genica | gnomAD | pLI, LOEUF (PVS1 en genes intolerantes a LoF) |
| Significancia clinica | ClinVar | veredicto + estrellas de oro |
| **IA missense** | **ESM-2** | **LLR evolutivo -> PP3/BP4** |
| **IA missense** | **AlphaMissense** | **patogenicidad 0-1 -> PP3/BP4** |

## Efecto demostrado

Con un perfil de fibrosis quistica extraido de la nota clinica (10 terminos HPO), **CFTR
sube del puesto #3-4 (solo variante) al #1**; BRCA1 -mayor score de variante pero sin
relacion fenotipica- baja al #2. El diseno multicapa evita falsos positivos: TP53 P72R y
APOE R176C salen "deletereas" por ESM-2 (LLR -3.3 y -6.5), pero su alta frecuencia poblacional
(regla BA1) y la ausencia de coincidencia fenotipica las mantienen correctamente como benignas
-ESM-2 no atropella la evidencia clinica/poblacional. Los cambios de aminoacido se anotan por
region+alelo del VCF para que el residuo mutante puntuado por ESM-2 coincida con el alelo real
del paciente (p. ej. BRCA1 A566E, EGFR L858R, TP53 P72R).

## Caso real: epilepsia infantil (GRCh37)

Sobre un VCF clínico real (varón, 10 meses, epilepsia; GRCh37/hg19; 1730 variantes → 871 PASS
únicas → 42 raras/impactantes → 16 en genes de epilepsia), la hipótesis principal fue **SCN1A
p.Cys1376Arg** (síndrome de Dravet). Aquí las capas de IA se complementan: **ESM-2 8M la marcó
tolerada por error (LLR +0.08), pero AlphaMissense la clasificó correctamente como probablemente
patogénica (0.999)**, en concordancia con SIFT/PolyPhen, la ausencia en gnomAD y el criterio PM5
(el codón Cys1376 ya alberga en ClinVar la variante C1376Y *patogénica* y C1376S *probablemente
patogénica* — verificado vía ClinVar E-utilities; el cambio del paciente, C1376R, es distinto y aún no reportado). Es la ilustración de por qué no se
usa un único predictor: AlphaMissense corrige el punto débil de ESM-2 8M en casos límite.

## Uso completo

```bash
python raredx_pipeline.py input.vcf \
       --assembly GRCh37 \              # build del VCF (GRCh38 por defecto)
       --clinical-note nota.txt \       # IA-1: HPO desde texto (o --hpo "HP:...")
       --esm \                          # IA-2: ESM-2 en missense
       --alphamissense \                # IA-3: AlphaMissense en missense (sin GPU)
       --out-prefix salida/paciente \
       --email tu@institucion.org
```

**Dependencias:** `requests` (base); `torch fair-esm` (para `--esm`); `anthropic` +
`ANTHROPIC_API_KEY` (para `--clinical-note`). `--alphamissense` **no requiere dependencias
extra ni GPU** (usa scores precalculados vía Ensembl VEP). Sin GPU, ESM-2 usa el modelo 8M en CPU.

## Salidas
- `<prefix>_report.html` - informe clinico: nota->HPO, ranking con LLR de ESM-2, fichas de evidencia
- `<prefix>_annotated.csv` - tabla completa con scores de variante, ESM-2, fenotipo y combinado

## Nota importante
Sistema de **apoyo a la decision, no un diagnostico**. Todo hallazgo debe confirmarse por
metodo ortogonal e interpretarse por un genetista clinico en el contexto completo del paciente.
Los scores de ESM-2 con el modelo 8M son una prueba de concepto; para uso serio, emplear
ESM-2 650M/3B (GPU) o AlphaMissense.
