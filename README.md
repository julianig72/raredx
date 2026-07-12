# raredx — Análisis de VCF con IA para diagnóstico de enfermedades raras

Pipeline moderno que anota un VCF y prioriza variantes candidatas para diagnóstico de
enfermedades raras, combinando **cuatro capas de evidencia clásica + cuatro capas de IA**, con
priorización dirigida por fenotipo al estilo de Exomiser. Soporta **GRCh38 y GRCh37/hg19**
(este último con liftover automático a GRCh38 donde hace falta).

Todo el análisis se apoya en **APIs REST públicas en vivo** (Ensembl, gnomAD, ClinVar, Open
Targets, OLS4): no hay que descargar bases de datos de gigabytes ni mantener índices locales.

## Documentación

Documentación científica detallada del pipeline (fuentes de datos, reglas ACMG con umbrales
exactos, predictores de IA, fórmulas de puntuación y capa agéntica):

- 📄 [Documentación científica (español)](docs/raredx_documentacion.html)
- 📄 [Scientific documentation (English)](docs/documentation_en.html)

Guía de despliegue de la herramienta web: [`web/README.md`](web/README.md).

## Instalación

Requiere **Python 3.11+**. GitHub Copilot es el proveedor LLM predeterminado para la nota clínica
y la capa agéntica; Anthropic queda disponible como fallback opcional.

```bash
git clone https://github.com/julianig72/raredx.git
cd raredx

# 1) Núcleo + proveedor GitHub Copilot
pip install -r requirements.txt
python -m copilot download-runtime
gh auth login   # cuenta con acceso a GitHub Copilot

# 2) Opcional — ESM-2 (IA-2). Añade ~2 GB (PyTorch). Sin GPU usa el modelo 8M en CPU.
pip install torch fair-esm

# 3) Opcional — fallback Anthropic
pip install anthropic
export ANTHROPIC_API_KEY="sk-ant-..."
```

| Necesitas… | Instala | Requiere clave/GPU |
|------------|---------|--------------------|
| Anotación + fenotipo + ranking | `requests` | no |
| AlphaMissense (IA-3) | *(nada extra)* — vía Ensembl VEP | no |
| ESM-2 (IA-2) | `torch fair-esm` | GPU opcional (8M corre en CPU) |
| HPO desde texto (IA-1) | `github-copilot-sdk` | sesión GitHub con Copilot |
| Capa agéntica (IA-4) | `github-copilot-sdk` | sesión GitHub con Copilot |
| Fallback LLM | `anthropic` | `ANTHROPIC_API_KEY` |

> **Conexión a internet obligatoria:** todas las anotaciones se resuelven contra APIs REST en vivo.

Para la **herramienta web** (interfaz para clínicos), ver [`web/README.md`](web/README.md).

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
                        Ranking de variantes
                                        v
        [IA-4] CAPA AGÉNTICA (estilo DeepRare, opcional)
        autorreflexión sobre candidatos (fenotipo + herencia)
        -> diferencial de ENFERMEDADES + enlaces verificados
                                        v
                        Ranking + informe HTML + CSV
```

## Las cuatro capas de IA

### IA-1 - Extraccion de fenotipo HPO desde nota clinica (LLM)
En vez de introducir codigos HPO a mano, **pegas la nota clinica en lenguaje natural** y un
LLM extrae los fenotipos presentes, que OLS4 aterriza a terminos HPO oficiales. Trabaja en
varios idiomas y captura descripciones coloquiales de signos, no solo terminologia medica
formal (p. ej. una descripcion de sudor salado se normaliza a `HP:0012236` "cloruro en sudor
elevado"; "dedos en palillo de tambor" a `HP:0100759`). Cada termino extraido queda ligado a
la frase de la nota que lo justifica, para trazabilidad.

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

### IA-4 - Capa agéntica: autorreflexión + diagnóstico diferencial trazable
Inspirada en **DeepRare** ("An agentic system for rare disease diagnosis with traceable
reasoning", *Nature* 2025, [doi:10.1038/s41586-025-10097-9](https://doi.org/10.1038/s41586-025-10097-9)). En lugar de entregar solo una lista rankeada de variantes,
un LLM **razona sobre su propia salida** en tres pasos:

1. **Autorreflexión.** Revisa cada candidato del top-K frente al fenotipo del paciente y al
   **patrón de herencia vs. la cigosidad observada**, emitiendo un veredicto `support`/`uncertain`/`refute`
   con su razonamiento. Detecta incoherencias que un ranking numérico no ve — p. ej. *una variante
   heterocigota única en un gen recesivo sin segundo golpe es un portador, no la causa*, o *una
   llamada heterocigota ligada al X en un varón es sospechosa*. Si **todos** los candidatos se
   descartan, amplía la ventana de búsqueda (K += paso) y reitera — el análogo fiel al bucle de
   DeepRare que aumenta la profundidad N y vuelve a iterar.
2. **Síntesis del diferencial.** Agrupa las variantes que sobreviven en **hipótesis a nivel de
   enfermedad** (enfermedad, gen(es), variante(s), herencia, razonamiento) y las ordena por
   verosimilitud.
3. **Verificación de enlaces (anti-alucinación).** Cada cita usa solo URLs **pre-construidas de
   forma determinista** (ClinVar, OMIM, Open Targets, Ensembl, PubMed) y se **comprueba que cada
   enlace resuelve** antes de incluirlo — el mismo mecanismo de trazabilidad de DeepRare. El LLM
   nunca inventa enlaces.

```bash
python raredx_pipeline.py input.vcf --agentic --alphamissense --assembly GRCh37 \
       --hpo "HP:0001250,HP:0011097" --out-prefix out/p
```

El resultado es que la capa **reordena por razonamiento clínico** en vez de solo por score
numérico: una variante puede tener el score más alto de la lista y aun así ser descartada por
incoherencia con la herencia (p. ej. un único alelo heterocigoto en un gen recesivo, o una
llamada ligada al X incompatible con el sexo), mientras que otra con menor score crudo pero
coherente con el fenotipo y la herencia asciende. El alcance es configurable con `start_k`
(por defecto 8): solo el top-K pasa por la autorreflexión, y la ventana se amplía únicamente
si todos los candidatos se refutan. El CSV marca cada fila con `agentic_evaluated` (yes/no).

## Motor de variante (ACMG-lite)

| Capa | Fuente | Aporta |
|------|--------|--------|
| Consecuencia funcional | Ensembl VEP | missense/stop/frameshift + SIFT/PolyPhen |
| Frecuencia poblacional | gnomAD r4 | AF (BA1/BS1: comun => benigna) |
| Restriccion genica | gnomAD | pLI, LOEUF como contexto; no bastan para activar PVS1 |
| Significancia clinica | ClinVar | veredicto + estrellas de oro |
| **IA missense** | **ESM-2** | **LLR evolutivo -> PP3/BP4** |
| **IA missense** | **AlphaMissense** | **patogenicidad 0-1 -> PP3/BP4** |

## Cómo se combina la evidencia (scoring)

El motor calcula un `variant_score` absoluto (0-1) por variante sumando pesos ACMG-lite; no lo
normaliza contra las demás variantes del mismo VCF. Las capas de IA missense participan antes de
la clasificación final:

- Una consecuencia `frameshift`, `stop_gained` o de splicing se marca como `LoF_predicted`, pero
  **no recibe PVS1 automáticamente**. PVS1 exige demostrar que la pérdida de función es el mecanismo
  de la enfermedad y revisar transcrito, NMD y región funcional; pLI/LOEUF por sí solos no bastan.
- La clase de cinco niveles usa combinaciones de fuerzas ACMG/AMP, no una suma de puntos. Incluso
  cuando PVS1 se valida independientemente, `PVS1 + PM2` corresponde a **Likely pathogenic**, no a
  **Pathogenic**.
- **ESM-2** y **AlphaMissense** aportan `PP3` (deletérea) o `BP4` (tolerada) y ajustan el
  `variant_score`; AlphaMissense pesa algo más por ser un predictor clínicamente calibrado. Para
  evitar contar varias veces evidencia computacional correlacionada, todas las etiquetas de la
  familia PP3 (o BP4) cuentan una sola vez en la clasificación.
- La **frecuencia poblacional manda sobre el in-silico**: una variante que un predictor marque
  "deletérea" pero que sea común en gnomAD (regla `BA1`/`BS1`) se mantiene benigna — la IA no
  atropella la evidencia poblacional ni la clínica de ClinVar.
- Los cambios de aminoácido se anotan **por región + alelo** del propio VCF (no por rsID), para
  que el residuo mutante que se puntúa corresponda exactamente al alelo del paciente y no a otro
  alelo del mismo locus.

Si se aporta fenotipo, el score final combina variante y fenotipo:

```
combined = 0.55 * variant_score + 0.45 * pheno_score     (× 0.5 si la variante no pasa el QC del VCF)
```

Así, entre dos variantes con evidencia molecular parecida, la que encaja con la clínica del
paciente asciende en el ranking — el principio de priorización dirigida por fenotipo de Exomiser.
El `pheno_score` se calcula para cada enfermedad por separado, multiplicando la cobertura HPO por
la fuerza de asociación de Open Targets, y se conserva la mejor enfermedad. Nunca se suman
fenotipos de enfermedades distintas del mismo gen. Un perfil con una sola HPO genera una
advertencia de baja especificidad.

## Uso completo

```bash
python raredx_pipeline.py input.vcf \
       --assembly GRCh37 \              # override opcional; por defecto se detecta desde el VCF
       --clinical-note nota.txt \       # IA-1: HPO desde texto (o --hpo "HP:...")
       --esm \                          # IA-2: ESM-2 en missense
       --alphamissense \                # IA-3: AlphaMissense en missense (sin GPU)
       --agentic \                      # IA-4: autorreflexión + diferencial trazable (LLM)
       --reflect-k 12 \                 # nº de candidatos que revisa la IA (por defecto 8)
       --father padre.vcf \             # trío: detecta de novo → criterio ACMG PS2
       --mother madre.vcf \             # trío: VCF de la madre
       --out-prefix salida/paciente \
       --email tu@institucion.org
```

**Análisis de trío (opcional).** Si se aportan los VCF de los progenitores (`--father`/`--mother`),
cada variante candidata se clasifica por herencia buscando el mismo alelo en los padres:
`de_novo` (el hijo porta el alelo y ambos progenitores tienen una llamada explícita sin él —por
ejemplo `0/0`— → añade el criterio ACMG **PS2**),
`paterna`/`materna`, `biparental` u `homozygous_recessive`. Las columnas `inheritance`,
`father_gt` y `mother_gt` se añaden al CSV, y el informe muestra una columna "Herencia".
La mera ausencia del sitio en un VCF parental de solo variantes no se considera prueba de
ausencia del alelo.

En un VCF con varias muestras es obligatorio indicar `--sample` con un identificador presente en
la cabecera; el pipeline nunca selecciona silenciosamente la primera muestra.

En la interfaz web, la expansión ontológica se muestra antes del análisis como una lista separada de
términos directos y ancestros. El usuario puede eliminar o añadir términos y el motor conserva por
separado los fenotipos directos para no convertir automáticamente todos los ancestros en coincidencias
directas. Durante el análisis se muestran tiempo transcurrido, ritmo, ETA y capas activas, y puede
solicitarse una cancelación segura desde la propia barra de progreso.

**Dependencias:** `requests` (base); `github-copilot-sdk` + una sesión GitHub autenticada con
Copilot (nota clínica y capa agéntica); `torch fair-esm` (para `--esm`). Anthropic puede configurarse
como fallback con `ANTHROPIC_API_KEY`. `--alphamissense` **no requiere dependencias extra ni GPU**
(usa scores precalculados vía Ensembl VEP). Sin GPU, ESM-2 usa el modelo 8M en CPU. Si no hay ningún
LLM disponible, el análisis se completa y muestra una advertencia.

## Salidas
- `<prefix>_report.html` - informe clinico: nota->HPO, ranking con LLR de ESM-2, fichas de evidencia
- `<prefix>_annotated.csv` - tabla completa con scores de variante, ESM-2, fenotipo y combinado

## Nota importante
Sistema de **apoyo a la decision, no un diagnostico**. Todo hallazgo debe confirmarse por
metodo ortogonal e interpretarse por un genetista clinico en el contexto completo del paciente.
Los scores de ESM-2 con el modelo 8M son una prueba de concepto; para uso serio, emplear
ESM-2 650M/3B (GPU) o AlphaMissense.
