# raredx web — herramienta clínica de análisis de VCF

Interfaz web para que un clínico analice un VCF sin tocar la línea de comandos: sube el
archivo, detecta automáticamente el ensamblaje, permite extraer HPO de una nota clínica y
revisarlos antes del análisis, y genera un informe priorizado descargable (HTML + CSV).

## Arquitectura

```
  Navegador (web/static/index.html)
     |  POST /api/extract-hpo  (nota clínica → HPO revisables)
     |  POST /api/expand-hpo   (perfil directo → ancestros HPO editables)
     |  multipart POST /api/analyze  (VCF + perfil HPO aprobado + opciones)
     v
  FastAPI (web/server.py)
     |  hilo en segundo plano -> raredx_pipeline.run_pipeline(...)
     |  progreso/ETA en vivo -> GET /api/status/{job}   (polling)
     |  cancelación cooperativa -> POST /api/cancel/{job}
     v
  Resultados por trabajo: <job>/result_report.html + result_annotated.csv
     |  GET /api/report/{job}   GET /api/csv/{job}
```

El servidor **reutiliza el mismo motor** que la CLI (`raredx_pipeline.run_pipeline`), de modo
que la web y la línea de comandos producen resultados idénticos.

## Ejecutar en local

```bash
pip install -r web/requirements.txt
python -m copilot download-runtime
# autenticar una cuenta con acceso a GitHub Copilot: gh auth login
# desde la raíz del repo:
uvicorn web.server:app --host 127.0.0.1 --port 8000
# abrir http://127.0.0.1:8000
```

## Ejecutar con Docker

```bash
docker build -f web/Dockerfile -t raredx-web .
docker run -p 8000:8000 -v raredx_data:/data \
       -e GH_TOKEN=token_de_usuario_con_copilot \
       -e RAREDX_EMAIL=tu@institucion.org \
       raredx-web
```

## Endpoints de la API

| Método | Ruta | Descripción |
|--------|------|-------------|
| `POST` | `/api/extract-hpo` | Extrae términos HPO revisables desde `clinical_note` y devuelve el proveedor LLM usado. |
| `POST` | `/api/expand-hpo` | Resuelve términos directos y devuelve sus ancestros HPO etiquetados para revisión. |
| `POST` | `/api/analyze` | Sube VCF + opciones; devuelve `{job_id, assembly}`. `assembly=auto` detecta GRCh37/38 desde el VCF. |
| `POST` | `/api/cancel/{job}` | Solicita la cancelación segura de un análisis activo. |
| `GET`  | `/api/status/{job}` | Estado, porcentaje, tiempo, ritmo, ETA y capas activas; al terminar incluye el top-15. |
| `GET`  | `/api/report/{job}` | Informe HTML completo. |
| `GET`  | `/api/csv/{job}` | Tabla anotada completa (CSV). |

## Variables de entorno

| Variable | Por defecto | Uso |
|----------|-------------|-----|
| `RAREDX_DATA_DIR` | temp del SO | Directorio de trabajos (montar volumen con política de borrado). |
| `RAREDX_MAX_MB` | `50` | Tamaño máximo de subida. |
| `RAREDX_MAX_ACTIVE_JOBS` | `4` | Máximo de análisis simultáneos o en cola. |
| `RAREDX_JOB_TTL_HOURS` | `24` | Retención antes de borrar trabajos terminados y directorios huérfanos. |
| `RAREDX_EMAIL` | — | Email de contacto para NCBI E-utilities (opcional). |
| `RAREDX_LLM_PROVIDER` | `auto` | `auto`, `copilot`, `anthropic` o `host`. |
| `RAREDX_COPILOT_MODEL` | `gpt-5-mini` | Modelo disponible en la suscripción de Copilot. |
| `RAREDX_LLM_TIMEOUT_SECONDS` | `120` | Timeout de cada inferencia LLM. |
| `RAREDX_MAX_NOTE_CHARS` | `20000` | Longitud máxima de una nota clínica. |
| `RAREDX_MAX_HPO_CHARS` | `10000` | Longitud máxima de la lista HPO enviada al análisis. |
| `RAREDX_MAX_SAMPLE_CHARS` | `200` | Longitud máxima del identificador de muestra. |
| `RAREDX_MAX_LLM_REQUESTS` | `2` | Extracciones HPO simultáneas. |
| `RAREDX_MAX_VARIANTS` | `50000` | Máximo de alelos llamados en la muestra que procesa un trabajo. |
| `RAREDX_JOB_TIMEOUT_SECONDS` | `3600` | Límite total de ejecución por análisis. |
| `GH_TOKEN` / `COPILOT_GITHUB_TOKEN` | sesión local | Autenticación GitHub para Copilot en servidores/contenedores. |
| `PORT` | `8000` | Puerto del servidor. |
| `ANTHROPIC_API_KEY` | — | Fallback opcional si Copilot no está disponible. |

## Antes de usar con pacientes reales — IMPORTANTE

- **Sin autenticación:** el servidor no la incluye. Colócalo detrás de un proxy inverso con
  autenticación (OAuth2/OIDC) y **TLS**. No expongas el puerto directamente.
- **Datos del paciente (PHI):** solo se guardan en el directorio del trabajo y se eliminan al
  vencer `RAREDX_JOB_TTL_HOURS` cuando el servicio procesa una nueva solicitud. Configura
  `RAREDX_DATA_DIR` en un volumen con política de retención/borrado acorde a tu normativa
  (RGPD / HIPAA según jurisdicción).
- **Proveedor LLM:** al usar nota clínica o razonamiento agéntico, el texto clínico y la evidencia
  seleccionada se envían al proveedor configurado (GitHub Copilot o Anthropic). Verifica que el uso
  y la retención del proveedor cumplan la política de datos clínicos de tu organización.
- **Apoyo a la decisión, no diagnóstico:** cada informe lo indica. Todo hallazgo debe
  confirmarse por método ortogonal (Sanger) e interpretarse por un genetista clínico.
- El VCF se **analiza, nunca se ejecuta**; las subidas están limitadas por tamaño.

## Pruebas

```bash
pip install httpx   # TestClient
python -m pytest web/test_server.py -v    # (o el script de humo incluido)
```
