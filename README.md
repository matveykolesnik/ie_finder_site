# ie_finder_site

Тонкий форк для сайта. Поиск элемента не скопирован: его выполняет `MGE_finder` ([rljech13/MGE_finder](https://github.com/rljech13/MGE_finder.git), коммит в `UPSTREAM`). Здесь только выкладка attL/attR на исходную сборку.

Рядом с этим каталогом должен лежать чекаут основной версии, либо путь задаётся явно:

```bash
./run.sh strain.fasta outdir
MGE_FINDER=/path/to/MGE_finder ./run.sh genomes_dir outdir
```

В `outdir/` три файла на сборку:

- `strain.ie.gff3`
- `strain.ie.gbk`
- `strain.ie.report.txt` — один отчёт: сколько интеграз и пар с тРНК нашлось, почему кандидат отброшен, координаты опубликованных элементов и непустые логи шагов

Координаты 1-based, на контигах входного FASTA: span элемента, attL, attR и CDS интегразы. Промежуточные таблицы и вырезанные острова пишутся во временный каталог и удаляются. `KEEP_WORK=1` оставляет его. `ANNOTATE_ALL=1` добавляет кандидатов с attL, которые не прошли более поздний фильтр. Если прогон упал, отчёт всё равно пишется из того, что успело посчитаться.

Пороги читаются из `MGE_finder/finder_pipeline/ie_finder_config.yaml`. Окружение — `MGE_finder/envs/IE_finder.yaml` (`USE_CONDA=1` на чистой машине).
