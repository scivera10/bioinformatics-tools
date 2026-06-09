#!/usr/bin/env bash
# =============================================================================
# extraccion_region_interes.sh
# Extrae lecturas de una región genómica de interés para todas las muestras,
# genera archivos FASTQ por par y guarda los resultados en directorios propios.
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# CONFIGURACIÓN — edita esta sección antes de ejecutar
# -----------------------------------------------------------------------------
DIRECTORIO_BASE="$(pwd)"

# Región genómica de interés en formato samtools: chr:inicio-fin
REGION="chrX:1000000-1005000"

# Lista de nombres de muestra
MUESTRAS=(
    Sample1
    Sample2
    Sample3
)

# Estructura esperada del BAM de entrada por muestra:
#   ${DIRECTORIO_BASE}/<MUESTRA>/resecuenciacion/2_bowtie2_mapping/<MUESTRA>.sorted.bam
# Directorio de salida por muestra:
#   ${DIRECTORIO_BASE}/NGS_SANGER/<MUESTRA>_SANGER/
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
    BAM_ENTRADA="${DIRECTORIO_BASE}/${MUESTRA}/resecuenciacion/2_bowtie2_mapping/${MUESTRA}.sorted.bam"

    # --- Directorio de salida exclusivo para esta muestra ---
    DIR_SALIDA="${DIRECTORIO_BASE}/NGS_SANGER/${MUESTRA}_SANGER"

    # --- Verificar que el BAM de entrada existe ---
    if [[ ! -f "${BAM_ENTRADA}" ]]; then
        echo "  [AVISO] No se encontró el BAM para ${MUESTRA}: ${BAM_ENTRADA}"
        echo "  Saltando muestra..."
        continue
    fi

    # --- Verificar que el directorio de salida no existe ya ---
    if [[ -d "${DIR_SALIDA}" ]]; then
        echo "  [AVISO] El directorio de salida ya existe: ${DIR_SALIDA}"
        echo "  Saltando muestra para evitar sobreescritura..."
        continue
    fi

    mkdir -p "${DIR_SALIDA}"

    # --- Rutas de archivos de salida ---
    BAM_REGION="${DIR_SALIDA}/region_interes.bam"
    BAM_NAMESORTED="${DIR_SALIDA}/region_interes.name_sorted.bam"
    FASTQ_R1="${DIR_SALIDA}/${MUESTRA}_R1.fq.gz"
    FASTQ_R2="${DIR_SALIDA}/${MUESTRA}_R2.fq.gz"
    FASTQ_SINGLETONS="${DIR_SALIDA}/singletons.fastq.gz"

    # --- Paso 1: Extraer la región de interés ---
    echo "  [1/5] Extrayendo región ${REGION}..."
    samtools view -b \
        "${BAM_ENTRADA}" \
        "${REGION}" \
        > "${BAM_REGION}"

    # --- Paso 2: Contar alineamientos extraídos ---
    echo "  [2/5] Contando alineamientos..."
    N_READS=$(samtools view -c "${BAM_REGION}")
    echo "        Alineamientos extraídos: ${N_READS}"

    # --- Paso 3: Estadísticas del BAM regional ---
    echo "  [3/5] Calculando flagstat..."
    samtools flagstat "${BAM_REGION}" \
        > "${DIR_SALIDA}/flagstat.txt"
    cat "${DIR_SALIDA}/flagstat.txt"

    # --- Paso 4: Ordenar por nombre de lectura ---
    echo "  [4/5] Ordenando por nombre de lectura..."
    samtools sort -n \
        -o "${BAM_NAMESORTED}" \
        "${BAM_REGION}"

    # --- Paso 5: Convertir a FASTQ ---
    echo "  [5/5] Convirtiendo a FASTQ..."
    samtools fastq \
        -1 "${FASTQ_R1}" \
        -2 "${FASTQ_R2}" \
        -s "${FASTQ_SINGLETONS}" \
        -0 /dev/null \
        -n \
        "${BAM_NAMESORTED}"

    echo "  [OK] ${MUESTRA} completada → ${DIR_SALIDA}"

done

echo ""
echo "======================================================"
echo "  Proceso finalizado para todas las muestras."
echo "======================================================"
