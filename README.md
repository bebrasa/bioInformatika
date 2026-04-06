# Домашнее задание 1 (Bioinformatics)

## 1) Фенотип и описание

**Фенотип:** наследственная предрасположенность к раку молочной железы и яичников  
**OMIM:** [Breast-ovarian cancer, familial, susceptibility](https://omim.org)

Кратко: это наследственный синдром, при котором из-за патогенных вариантов в генах репарации ДНК существенно повышается риск раннего развития рака молочной железы и/или яичников. Наследование чаще аутосомно-доминантное по предрасположенности, а молекулярной основой часто служит нарушение пути гомологичной рекомбинации.

## 2) Ассоциированные гены

- `BRCA1`
- `BRCA2`

## 3) FASTA-файлы с последовательностями генов (NCBI Nucleotide)

Сравнивались человек (`Homo sapiens`) и модельный организм мышь (`Mus musculus`).

- `data/fasta/BRCA1_Homo_sapiens_NM_007294.4.fasta`
- `data/fasta/Brca1_Mus_musculus_NM_009764.fasta`
- `data/fasta/BRCA2_Homo_sapiens_NM_000059.4.fasta`
- `data/fasta/Brca2_Mus_musculus_NM_009765.fasta`

## 4) Файлы с парными выравниваниями

Использованы два инструмента/подхода:
- Needleman-Wunsch (глобальное выравнивание)
- Smith-Waterman (локальное выравнивание)

Результаты:
- `results/alignments/BRCA1_global_needleman_wunsch.txt`
- `results/alignments/BRCA1_local_smith_waterman.txt`
- `results/alignments/BRCA2_global_needleman_wunsch.txt`
- `results/alignments/BRCA2_local_smith_waterman.txt`

Скрипт для воспроизведения:
- `scripts/run_alignments.py`

Запуск:

```bash
.venv/bin/python scripts/run_alignments.py
```

## 5) Оценка качества выравниваний и выбор лучшего

### BRCA1 (человек vs мышь)
- **Global (Needleman-Wunsch):** identity = **75.93%**, aligned non-gap = 6511, score = 7487.0
- **Local (Smith-Waterman):** identity = **75.81%**, aligned non-gap = 6507, score = 7520.0

### BRCA2 (человек vs мышь)
- **Global (Needleman-Wunsch):** identity = **75.65%**, aligned non-gap = 10651, score = 11842.0
- **Local (Smith-Waterman):** identity = **75.41%**, aligned non-gap = 10621, score = 12041.0

### Вывод

Для данной задачи (сравнение **ортологических генов целиком**) более информативным считаю **глобальное выравнивание Needleman-Wunsch**:
- оно оптимизирует соответствие по всей длине последовательностей;
- даёт немного более высокую долю идентичности для обоих генов;
- лучше подходит для итогового межвидового сравнения полноразмерных транскриптов.

Локальное выравнивание Smith-Waterman полезно как дополнительный контроль, чтобы выделять наиболее консервативные участки, но как основное для сравнения целых последовательностей уступает глобальному подходу.
