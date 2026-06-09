#!/usr/bin/env bash
# =============================================================================
# alineamiento_sanger.sh
# Alinea las lecturas de la región de interés de cada muestra contra la
# secuencia consenso Sanger usando Bowtie2. Genera un BAM ordenado e indexado
# y un log de alineamiento por muestra.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# CONFIGURACIÓN — edita esta sección antes de ejecutar
# -----------------------------------------------------------------------------
DIRECTORIO_BASE="."                        # Directorio raíz del proyecto
DIR_REFERENCIA="${DIRECTORIO_BASE}/sanger_referencia"
INDICE_SANGER="${DIR_REFERENCIA}/sanger_index"   # Prefijo del índice Bowtie2 (sin extensión)
THREADS=4                                  # Núcleos para Bowtie2

# Lista de nombres de muestra (una por línea, sin espacios)
MUESTRAS=(
    Sample1
    Sample2
    Sample3
)

# Estructura esperada de directorios por muestra:
#   ${DIRECTORIO_BASE}/<MUESTRA>_SANGER/<MUESTRA>_R1.fq.gz
#   ${DIRECTORIO_BASE}/<MUESTRA>_SANGER/<MUESTRA>_R2.fq.gz
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# BUCLE PRINCIPAL
# -----------------------------------------------------------------------------
for MUESTRA in "${MUESTRAS[@]}"; do

    echo ""
    echo "======================================================"
    echo "  Procesando muestra: ${MUESTRA}"
    echo "======================================================"

    # --- Rutas de entrada ---
    DIR_MUESTRA="${DIRECTORIO_BASE}/${MUESTRA}_SANGER"
    FASTQ_R1="${DIR_MUESTRA}/${MUESTRA}_R1.fq.gz"
    FASTQ_R2="${DIR_MUESTRA}/${MUESTRA}_R2.fq.gz"

    # --- Rutas de salida ---
    BAM_CRUDO="${DIR_MUESTRA}/${MUESTRA}.sanger.bam"
    BAM_SORTED="${DIR_MUESTRA}/${MUESTRA}.sanger.sorted.bam"
    LOG_BOWTIE2="${DIR_MUESTRA}/${MUESTRA}.sanger.bowtie2.log"

    # --- Verificar que los FASTQ existen ---
    if [[ ! -f "${FASTQ_R1}" || ! -f "${FASTQ_R2}" ]]; then
        echo "  [AVISO] No se encontraron los FASTQ para ${MUESTRA}."
        echo "  Saltando muestra..."
        continue
    fi

    # --- Verificar que el BAM final no existe ya ---
    if [[ -f "${BAM_SORTED}" ]]; then
        echo "  [AVISO] El BAM de salida ya existe: ${BAM_SORTED}"
        echo "  Saltando muestra para evitar sobreescritura..."
        continue
    fi

    # --- Paso 1: Alineamiento con Bowtie2 ---
    echo "  [1/3] Alineando contra la referencia Sanger..."
    bowtie2 \
        -x "${INDICE_SANGER}" \
        -1 "${FASTQ_R1}" \
        -2 "${FASTQ_R2}" \
        --no-unal \
        -p "${THREADS}" \
        2> "${LOG_BOWTIE2}" \
    | samtools view -b -o "${BAM_CRUDO}"

    cat "${LOG_BOWTIE2}"

    # --- Paso 2: Ordenar por coordenada ---
    echo "  [2/3] Ordenando BAM por coordenada..."
    samtools sort \
        -o "${BAM_SORTED}" \
        "${BAM_CRUDO}"

    # --- Paso 3: Indexar el BAM final ---
    echo "  [3/3] Indexando BAM..."
    samtools index "${BAM_SORTED}"

    # --- Eliminar BAM crudo intermedio ---
    rm "${BAM_CRUDO}"

    echo "  [OK] ${MUESTRA} completada → ${BAM_SORTED}"

done

echo ""
echo "======================================================"
echo "  Alineamiento contra Sanger finalizado."
echo "======================================================"
