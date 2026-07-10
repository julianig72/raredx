# raredx web — herramienta clínica de análisis de VCF

Interfaz web para que un clínico analice un VCF sin tocar la línea de comandos: sube el
archivo, opcionalmente pega una nota clínica y/o términos HPO, elige el ensamblaje y las
capas de IA, y obtiene un informe priorizado descargable (HTML + CSV).

## Arquitectura

```
  Navegador (web/static/index.html)
     |  multipart POST /api/analyze  (VCF + nota + HPO + opciones)
     v
  FastAPI (web/server.py)
     |  hilo en segundo plano -> raredx_pipeline.run_pipeline(...)
     |  progreso en vivo -> GET /api/status/{job}   (polling)
     v
  Resultados por trabajo: <job>/result_report.html + result_annotated.csv
     |  GET /api/report/{job}   GET /api/csv/{job}
```

El servidor **reutiliza el mismo motor** que la CLI (`raredx_pipeline.run_pipeline`), de modo
que la web y la línea de comandos producen resultados idénticos.

## Ejecutar en local

```bash
pip install -r web/requirements.txt
# desde la raíz del repo:
uvicorn web.server:app --host 127.0.0.1 --port 8000
# abrir http://127.0.0.1:8000
```

## Ejecutar con Docker

```bash
docker build -f web/Dockerfile -t raredx-web .
docker run -p 8000:8000 -v raredx_data:/data \
       -e RAREDX_EMAIL=tu@institucion.org \
       raredx-web
```

## Endpoints de la API

| Método | Ruta | Descripción |
|--------|------|-------------|
| `POST` | `/api/analyze` | Sube VCF + opciones; devuelve `{job_id}`. Campos: `vcf` (archivo), `sample`, `hpo`, `clinical_note`, `assembly` (GRCh38/GRCh37), `esm`, `alphamissense`. |
| `GET`  | `/api/status/{job}` | Estado y progreso; al terminar incluye el top-15 de candidatos. |
| `GET`  | `/api/report/{job}` | Informe HTML completo. |
| `GET`  | `/api/csv/{job}` | Tabla anotada completa (CSV). |

## Variables de entorno

| Variable | Por defecto | Uso |
|----------|-------------|-----|
| `RAREDX_DATA_DIR` | temp del SO | Directorio de trabajos (montar volumen con política de borrado). |
| `RAREDX_MAX_MB` | `50` | Tamaño máximo de subida. |
| `RAREDX_EMAIL` | — | Email de contacto para NCBI E-utilities (opcional). |
| `PORT` | `8000` | Puerto del servidor. |
| `ANTHROPIC_API_KEY` | — | Necesaria para la capa de nota clínica (LLM). |

## Antes de usar con pacientes reales — IMPORTANTE

- **Sin autenticación:** el servidor no la incluye. Colócalo detrás de un proxy inverso con
  autenticación (OAuth2/OIDC) y **TLS**. No expongas el puerto directamente.
- **Datos del paciente (PHI):** solo se guardan en el directorio del trabajo. Configura
  `RAREDX_DATA_DIR` en un volumen con política de retención/borrado acorde a tu normativa
  (RGPD / HIPAA según jurisdicción).
- **Apoyo a la decisión, no diagnóstico:** cada informe lo indica. Todo hallazgo debe
  confirmarse por método ortogonal (Sanger) e interpretarse por un genetista clínico.
- El VCF se **analiza, nunca se ejecuta**; las subidas están limitadas por tamaño.

## Pruebas

```bash
pip install httpx   # TestClient
python -m pytest web/test_server.py -v    # (o el script de humo incluido)
```
