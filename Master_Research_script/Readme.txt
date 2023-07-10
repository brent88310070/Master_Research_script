Data processing:
    資料處理可以分成三部分(SHARE, HCC, MB) 分別以1.1, 1.2, 1.3處理
    第三筆資料使用MIRA套件畫UMAP(1.3.1)
    第二、三筆資料需要多跑 ChromatinLinkGene.R (第一筆寫在同個檔案)

Infer cellular direction:
    Infer_cell_direction.R 推細胞群的分化方向
    所需檔案為從Data processing跑完，具有兩群細胞種類的*.rds檔

Find fate commitment genes:
    為研究延伸部分，分成兩個scripts
    1.FateDecisionGene.R 需要mouse_tf_targetGene.tsv資料
    2.ExpressionSmoothedByPseudotime
