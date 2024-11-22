#!/usr/bin/env nextflow

// get tumor and normal sample ids from bam header
process SampleID {
    container "${params.container_cgpmap}"

    input:
    path files

    output:
    env(normal_file), emit: normal_file
    env(normal_SM), emit: normal_SM
    env(tumor_file), emit: tumor_file
    env(tumor_SM), emit: tumor_SM

    script:
    """
    echo "listing items in files: ${files}"
    ls ${files}
    echo ""

    # file_1 RG
    rg_1=\$(samtools view -H ${files}/${params.file_1} | grep '^@RG' | head -1)
    
    # file_2 RG
    rg_2=\$(samtools view -H ${files}/${params.file_2} | grep '^@RG' | head -1)

    # get tumor normal designation for each file (bam)
    type_1=\$(samtools view -H ${files}/${params.file_1} | grep '^@RG' | head -1 | sed 's/.*DS:\\(.*\\)/\\1/' | awk -F'|' '{print \$NF}')
    type_2=\$(samtools view -H ${files}/${params.file_2} | grep '^@RG' | head -1 | sed 's/.*DS:\\(.*\\)/\\1/' | awk -F'|' '{print \$NF}')

    echo "rg_1: \${rg_1}"
    echo "type_1: \${type_1}"
    echo ""
    echo "rg_2: \${rg_2}"
    echo "type_2: \${type_2}"
    echo ""
    
    normal_file=""
    normal_SM=""
    tumor_file=""
    tumor_SM=""
    
    # assign proper files to 'normal_file and normal_SM' and 'tumor_file and tumor_SM'
    if [[ "\${type_1}" == "Normal" ]]; then
        normal_file=${params.file_1}
        normal_SM=\$(samtools view -H ${files}/${params.file_1} | grep '^@RG' | head -1 | grep -oP 'SM:\\K[^\\s]+')
    elif [[ "\${type_1}" == "Tumour" ]]; then
        tumor_file=${params.file_1}
        tumor_SM=\$(samtools view -H ${files}/${params.file_1} | grep '^@RG' | head -1 | grep -oP 'SM:\\K[^\\s]+')
    else
        echo "${files}/${params.file_1} type designation is neither 'Normal' nor 'Tumour'"
    fi

    if [[ "\${type_2}" == "Normal" ]]; then
        normal_file=${params.file_2}
        normal_SM=\$(samtools view -H ${files}/${params.file_2} | grep '^@RG' | head -1 | grep -oP 'SM:\\K[^\\s]+')
    elif [[ "\${type_2}" == "Tumour" ]]; then
        tumor_file=${params.file_2}
        tumor_SM=\$(samtools view -H ${files}/${params.file_2} | grep '^@RG' | head -1 | grep -oP 'SM:\\K[^\\s]+')
    else
        echo "${files}/${params.file_2} type designation is neither 'Normal' nor 'Tumour'"
    fi

    echo "normal_file: \${normal_file}"
    echo "normal_SM: \${normal_SM}"
    echo ""
    echo "tumor_file: \${tumor_file}"
    echo "tumor_SM: \${tumor_SM}"
    """
}