```bash
ls -1 /usr/share/data-minor-bioinf/assembly/* | xargs -tI{} ln -s {}
seqtk sample -s1124 oil_R1.fastq 5000000 > sub1.fastq
seqtk sample -s1124 oil_R2.fastq 5000000 > sub2.fastq
seqtk sample -s1124 oilMP_S4_L001_R1_001.fastq 1500000 > mp1.fastq
seqtk sample -s1124 oilMP_S4_L001_R2_001.fastq 1500000 > mp2.fastq
ls -1 | xargs -P 4 -tI{} fastqc {}
multiqc -o multiqc fastqc
platanus_trim sub1.fastq sub2.fastq
platanus_internal_trim mp1.fastq mp2.fastq
mkdir trimmed
ls -1 | xargs -P 4 -tI{} fastqc {}
cd report_trimmed
multiqc . -o .
platanus assemble -o Poil -t 4 -m 15 -f sub1.fastq.trimmed sub2.fastq.trimmed 2> assemble.log
platanus scaffold -o Poil -t 4 -c Poil_contig.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 mp1.fastq.int_trimmed mp2.fastq.int_trimmed 2> scaffold.log
platanus gap_close -o Poil -t 4 -c Poil_scaffold.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 mp1.fastq.int_trimmed mp2.fastq.int_trimmed 2>gapclose.log
echo scaffold1_cov231 > tmp.txt
seqtk subseq Poil_gapClosed.fa tmp.txt > oil_genome.fna
```

```bash
ls -1 /usr/share/data-minor-bioinf/assembly/* | xargs -tI{} ln -s {}
seqtk sample -s1124 oil_R1.fastq 4000000 > sub1.fastq
seqtk sample -s1124 oil_R2.fastq 4000000 > sub2.fastq
seqtk sample -s1124 oilMP_S4_L001_R1_001.fastq 1200000 > mp1.fastq
seqtk sample -s1124 oilMP_S4_L001_R2_001.fastq 1200000 > mp2.fastq
ls -1 | xargs -P 4 -tI{} fastqc {}
multiqc -o multiqc fastqc
platanus_trim sub1.fastq sub2.fastq
platanus_internal_trim mp1.fastq mp2.fastq
ls -1 | xargs -P 4 -tI{} fastqc {} -o ../report_trimmed
multiqc . -o .
platanus assemble -o Poil -t 4 -m 15 -f sub1.fastq.trimmed sub2.fastq.trimmed 2> assemble.log
platanus scaffold -o Poil -t 4 -c Poil_contig.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 mp1.fastq.int_trimmed mp2.fastq.int_trimmed 2> scaffold.log
platanus gap_close -o Poil -t 4 -c Poil_scaffold.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 mp1.fastq.int_trimmed mp2.fastq.int_trimmed 2>gapclose.log
echo scaffold1_cov184 > tmp.txt
seqtk subseq Poil_gapClosed.fa tmp.txt > oil_genome.fna
```

```python
import re
# Для читаемости создал класс в котором буду хранить скаффолды и контиги
class contig:
  def __init__(self, length, sequence):
    # Для задания и читаемости достаточно хранить длину и последовательность
    self.length = length
    self.sequence = sequence
    
def report(file): 
  # Создаем массив для хранения считанных контигов 
  contigs = []
  sum = 0

  with open(file, 'r') as f:
    for line in f:
      # Если строка начинается с > - значит перед нами новый континг. Создаем пустой объект и добавляем его в массив
      if line[0] == '>':
        contigs.append(contig(0, ""))
      else:
        # Обычную строчку сначала очищаем от спец. символов
        _tmp = line.replace('\r', '').replace('\n', '')
        # Прибавляем последовательность к последнему созданному контигу (можем так сделать из-за структуры файла)
        contigs[-1].sequence += _tmp
        # Прибавляем длину строки к длине последнего контига
        contigs[-1].length += len(_tmp)
        # Сразу считаем суммарную длину
        sum+=len(_tmp)

  # Сортируем по длине для рассчета N50
  contigs.sort(reverse=True, key = lambda x: x.length)

  N50 = 0
  # Идем по контигам
  for i in contigs:
    # Добавляем длину
    N50 += i.length
    # Если с добавлением длины текущего контига мы превышаем половину суммарной длины - печатаем статистику и записываем эту длину как N50
    if N50 > sum/2:
      print(f"File: " + file)
      print(f"Number of contigs: {len(contigs)}")
      print(f"Sum length: {sum}") 
      print(f"Maximum contig length: {contigs[0].length}")
      print(f"N50: {i.length}")
      # Возвращаем все считанные контиги, пригодится для поиска гэпов
      return contigs

def report_gaps(sequence):
  # Находим все подстроки содержащие N
  gaps = re.findall("N+", sequence)
  # Длина массива будет их количеством
  print(f"Number of gaps: {len(gaps)}")
  # Сумма длин - общей длиной пропусков
  print(f"Length of gaps: {sum([len(i) for i in gaps])}")
  # Возвращаем все найденные гэпы
  return gaps

contigs = report('m_Poil_contig.fa')
scaffolds = report('m_Poil_scaffold.fa')
gaps = report_gaps(scaffolds[0].sequence)
no_gaps = report('m_Poil_gapClosed.fa')
gaps = report_gaps(no_gaps[0].sequence)

# Префикс о означает optional
contigs_o = report('o_Poil_contig.fa')
scaffolds_o = report('o_Poil_scaffold.fa')
gaps_o = report_gaps(scaffolds_o[0].sequence)
no_gaps_o = report('o_Poil_gapClosed.fa')
gaps_o = report_gaps(no_gaps_o[0].sequence)
```
