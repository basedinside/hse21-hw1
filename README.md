# Обязательная часть
- С помощью команды seqtk выбираем случайно 5 миллионов чтений типа paired-end и 1.5 миллиона чтений типа mate-pairs (чтобы у каждого получился свой уникальный результат)
```bash
seqtk sample -s1124 oil_R1.fastq 5000000 > sub1.fastq
seqtk sample -s1124 oil_R2.fastq 5000000 > sub2.fastq
seqtk sample -s1124 oilMP_S4_L001_R1_001.fastq 1500000 > mp1.fastq
seqtk sample -s1124 oilMP_S4_L001_R2_001.fastq 1500000 > mp2.fastq
```
- С помощью программы fastQC и multiQC оценить качество исходных чтений и получить по ним общую статистику
```bash
ls -1 | xargs -P 4 -tI{} fastqc {}
multiqc fastqc
```
- С помощью программ platanus_trim и platanus_internal_trim подрезать чтения по качеству и удалить праймеры
```bash
platanus_trim sub1.fastq sub2.fastq
platanus_internal_trim mp1.fastq mp2.fastq
```
- С помощью программы fastQC и multiQC оценить качество подрезанных чтений и получить по ним общую статистику
```bash
ls -1 | xargs -P 4 -tI{} fastqc {}
multiqc fastqc
```
- С помощью программы “platanus assemble” собрать контиги из подрезанных чтений
```bash
platanus assemble -o Poil -t 4 -m 15 -f sub1.fastq.trimmed sub2.fastq.trimmed 2> assemble.log
```
- Написать код (в Jupiter ноутбуке и в Google Colab) для анализа полученных контигов (общее кол-во контигов, их общая длина, длина самого длинного контига, N50)
```python
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
```
![image](https://user-images.githubusercontent.com/71254839/139060801-257e1e0b-d563-4bef-bac4-d47435a52a6b.png)
- С помощью программы “platanus scaffold” собрать скаффолды из контигов, а также из подрезанных чтений
```bash
platanus scaffold -o Poil -t 4 -c Poil_contig.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 mp1.fastq.int_trimmed mp2.fastq.int_trimmed 2> scaffold.log
```
- Написать код (в Jupiter ноутбуке и в Google Colab) для анализа полученных скаффолдов (общее кол-во скаффолдов, их общая длина, длина самого длинного скаффолда, N50)
![image](https://user-images.githubusercontent.com/71254839/139063686-9b0f9c2c-18e2-4e82-8fba-f4072d0b6fe4.png)
- Для самого длинного скаффолда посчитать количество гэпов (участков, состоящих из букв NNNN) и их общую длину
```python
def report_gaps(sequence):
  # Находим все подстроки содержащие N
  gaps = re.findall("N+", sequence)
  # Длина массива будет их количеством
  print(f"Number of gaps: {len(gaps)}")
  # Сумма длин - общей длиной пропусков
  print(f"Length of gaps: {sum([len(i) for i in gaps])}")
  # Возвращаем все найденные гэпы
  return gaps
```
![image](https://user-images.githubusercontent.com/71254839/139063075-f4e48ed6-9aa8-4275-8dfe-a611c559d68f.png)
- С помощью программы “platanus gap_close” уменьшить кол-во гэпов с помощью подрезанных чтений
```bash
platanus gap_close -o Poil -t 4 -c Poil_scaffold.fa -IP1 sub1.fastq.trimmed sub2.fastq.trimmed -OP2 mp1.fastq.int_trimmed mp2.fastq.int_trimmed 2>gapclose.log
```
- Для самого длинного скаффолда посчитать количество гэпов (участков, состоящих из букв NNNN) и их общую длину
![image](https://user-images.githubusercontent.com/71254839/139063471-008a0f2b-8c9a-4a47-9890-5ce6025fe8e0.png)
# Дополнительная часть
- Попробовать собрать геном из меньшего кол-ва чтений и посмотреть как качество сборки (кол-во скаффолдов и гэпов) меняется с кол-вом исходных чтений
```bash
seqtk sample -s1124 oil_R1.fastq 4000000 > sub1.fastq
seqtk sample -s1124 oil_R2.fastq 4000000 > sub2.fastq
seqtk sample -s1124 oilMP_S4_L001_R1_001.fastq 1200000 > mp1.fastq
seqtk sample -s1124 oilMP_S4_L001_R2_001.fastq 1200000 > mp2.fastq
```
Я решил взять на 20% меньше, чем в основной части. Остальной код аналогичный основной части.
Я получил следующие данные:
- Для контигов
![image](https://user-images.githubusercontent.com/71254839/139064587-b1e87159-112e-4c65-aa16-92c5186507db.png)
- Для скаффолдов
![image](https://user-images.githubusercontent.com/71254839/139065365-595f4685-1bc7-41a8-b51d-91d56309046b.png)
- Для скаффолдов с закрытыми гэпами
![image](https://user-images.githubusercontent.com/71254839/139065328-7742e3fc-69fe-42d3-adcb-4559c5843f2c.png)

Можем заметить, что результаты получились неоднозначные.
