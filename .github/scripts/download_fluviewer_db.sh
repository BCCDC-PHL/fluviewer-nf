#!/bin/bash

mkdir -p .github/data/fluviewer_db

wget -O .github/data/fluviewer_db/FluViewer_db_v_0_1_8.fa.gz https://raw.githubusercontent.com/KevinKuchinski/FluViewer/main/FluViewer_db_v_0_1_8.fa.gz

gunzip .github/data/fluviewer_db/FluViewer_db_v_0_1_8.fa.gz
