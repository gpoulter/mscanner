"""Test data for genedrug_test.py"""

drugs_table="""PA10000\t17 beta-estradiol\t\t
PA10007\talbumin human\t\tAlbuminar-25|Albuminar-5|Albutein 25%|Albutein 5%|Buminate 25%|Buminate 5%|Plasbumin-25|Plasbumin-5|
PA10009\talefacept\t\tAmevive|
PA1001\t1-methyl-4-phenylpyridinium (MPP+)\t\t
PA10010\talemtuzumab\t\tCampath|
PA10011\talfacalcidol\t\tOne-Alpha|
PA10012\talteplase, recombinant\t\tActivase|Activase rt-PA|Cathflo Activase|"""

fulldrugs_correct={'PA10000': ['17 beta-estradiol'],
 'PA10007': ['albumin human',
             'Albuminar-25',
             'Albuminar-5',
             'Albutein 25%',
             'Albutein 5%',
             'Buminate 25%',
             'Buminate 5%',
             'Plasbumin-25',
             'Plasbumin-5'],
 'PA10009': ['alefacept', 'Amevive'],
 'PA1001': ['1-methyl-4-phenylpyridinium (MPP+)'],
 'PA10010': ['alemtuzumab', 'Campath'],
 'PA10011': ['alfacalcidol', 'One-Alpha'],
 'PA10012': ['alteplase, recombinant',
             'Activase',
             'Activase rt-PA',
             'Cathflo Activase']}

shortdrugs_correct={'PA10000': ['beta estradiol'],
 'PA10007': ['albumin human',
             'albuminar',
             'albuminar',
             'albutein',
             'albutein',
             'buminate',
             'buminate',
             'plasbumin',
             'plasbumin'],
 'PA10009': ['alefacept', 'amevive'],
 'PA1001': ['methyl phenylpyridinium mpp'],
 'PA10010': ['alemtuzumab', 'campath'],
 'PA10011': ['alfacalcidol', 'one alpha'],
 'PA10012': ['alteplase recombinant',
             'activase',
             'activase rt pa',
             'cathflo activase']}

genes_correct=[('HGF', 142, 145), ('Hepatocyte Growth\nFactorase', 113, 140), ('MAPK', 61, 65), ('mitogen-activated protein kinase', 27, 59)]

drugs_text=[
    "17 beta estradiol, Albuminar 25, Amevive",
    "methyl-4-phenylpyridinium (MPP+),alemtuzumab, Activate rt-PA",
    "one:Alpha, Activase-rt*PA"
    ]

drugs_correct=[
    {'PA10000': '17 beta-estradiol','PA10007': 'Albuminar-5','PA10009': 'Amevive'},
    {'PA1001': '1-methyl-4-phenylpyridinium (MPP+)', 'PA10010': 'alemtuzumab'},
    {'PA10012': 'Activase rt-PA', 'PA10011': 'One-Alpha'}
    ]

genefinder_cache={
    "'We observed an increase in mitogen-activated protein kinase (MAPK)\\nactivity when administering Amevive. We found Hepatocyte Growth\\nFactorase (HGF) activity decreased when administering Campath.'":
    [('HGF',142,145,1),
     ('Hepatocyte Growth\nFactorase',113,140,1),
     ('MAPK',61,65,1),
     ('mitogen-activated protein kinase',27,59,1)]}

genedrug_text="""We observed an increase in mitogen-activated protein kinase (MAPK)
activity when administering Amevive. We found Hepatocyte Growth
Factorase (HGF) activity decreased when administering Campath."""

genedrug_correct={
    'Amevive': set(['MAPK', 'mitogen-activated protein kinase']),
    'Campath': set(['HGF', 'Hepatocyte Growth\nFactorase'])
    }

test_articles=[
    Article(
    pmid=1,
    title="We observed an increase in mitogen-activated protein kinase (MAPK) activity when administering Amevive.",
    abstract="We found Hepatocyte Growth Factorase (HGF) activity decreased when administering Campath.",
    meshterms=["T1"]
    ),
    Article(
    pmid=2,
    title="Gobbeldy-gook",
    abstract="Gobbeldy-gook",
    meshterms=["T1"]
    )
    ]

test_cache={
    1:{'PA10010': set(['HGF', 'MAPK', 'Hepatocyte Growth Factorase', 'mitogen-activated protein kinase']),
                'PA10009': set(['HGF', 'MAPK', 'Hepatocyte Growth Factorase', 'mitogen-activated protein kinase'])},
    2: {}
    }

### THE FOLLOWING IS AN EXHAUSTIVE TEST OF THE DRUG LISTER

positive_articles="""11809184
12069159
9744524
1960624
8644731
10689270
9103127
11526473
11825235
9157990
9952058
11259360
2831655
10521338
9882456
9441946
10064574
11925851
9841604
12172398
9352571
9352573
11907494
11875365
10553725
10916115
7924124
10984540
9950666
12123487
9261824
9777427
10824074
12112247
9549641
10904024
2020918
8004130
9224780
10755579
9276194
2223426
10660000
8747411
9719084
10656877
9569045
9660842
3816023
11597824
9857976
7690693
10379516
2568832
8742240
2311333
11740344
8110777
12006904
2897240
12172218
1867960
11950793
9867757
9211656
11523725
8685072
10997944
10868744
10668853
11089838
10764140
8752164
9632445
10620209
12052139
9345314
9806540
11580286
9806549
8986921
12172211
12175731
7862671
12049175
6641087
12152005
11942593
3335642
7891353
12082588
12082589
9869506
8643688
11919084
2849209
9510461
11668216
9825838
9663807
10988267
9230185
12170773
9435432
11668219
12351588
2015731
10470083
8591723
10363930
8807659
10944550
10373227
8895751
12010835
7398195
11053256
2327982
9834040
12142730
12142731
9652563
11229339
10027427
12042667
11304783
9568693
11543872
10739171
8100166
10961881
12083975
11237013
9547362
11950782
12172213
8555090
10073976
11773861
9476037
8605025
9270093
9190321
12152008
10821369
7806690
12386645
12386646
12386647
10792337
10880412
11434505
9689128
9918138
9149697
10636265
10369259
9324115
10066160
12228189
10093988
11740338
10515402
11913724
11823760
11823761
9831561
8823236
12086014
11692081
8538349
9241661
11475331
9492384
11947894
9634537
3467886
11279528
8941024
10082212
9655391
9672899
7586933
12083394
9202749
7562475
9473242
7891318
11990384
8866821
11990382
11990381
11773859
11302935
12200397
11956513
9005998
3956053
9780130
11814344
12107857
12202283
6519160
1587064
11908756
9690230
11180734
11678778
12203989
9106238
9713905
10430787
11906176
10073515
10820139
12006396
9190854
9164836
8620056
9114722
10613621
10580024
12235454
12235455
9433390
9777180
2536896
8351162
11526468
1971708
10889553
9164413
12189368
11919305
8561894
10212557
9661139
9164419
11549275
11368430
10190398
4007030
9010701
10867985
7474246
11264452
8257179
9329515
9872059
7619675
11259359
8618881
12072547
9682271
10535720
8359181
7640149
1955133
3123997
12196913
71400
1857341
10898111
1785203
12004990
1675794
2195556
8521555
11159808
1782973
11446878
9386136
10903223
9551410
12145682
10736278
11875363
12108579
11875366
9802529
9269699
9731720
10704169
2060250
9484992
1505606
12172219
12172216
12172215
10801254
10895986
7491709
3598909
9726242
8143391
8952600
9407415
9613951
7605354
10790202
8772238
8094196
8242617
2104790
10647091
12042670
12042672
9629566
421578
12167489
12087351
12087352
12082591
11926893
10862524
9131488
8742444
9202739
9265947
10893495
11966400
9466980
11966404
11966407
10999650
10219239
11934710
9069453
12086290
12042669
10340917
9823305
9399966
7903454
10216075
9578183
9653159
10716719
12142724
1354642
9280407
9789606
12142728
11337938
9194910
11910656
11956512
8565783
9420339
10860023
11551222
2007317
10884039
11997281
12171903
11773863
12164772
12049181
11773860
11773867
11773866
10945315
11773864
9521254
11070098
9377473
2570698
10335726
9789061
11773862
9114822
11916544
9806904
10683861
7593601
11970993
8513845
9255563
8818577
11714868
9316860
9295054
12269967
8819499
9034160
10192756
9838066
11279519
12065207
9277043
10051705
10406363
9259381
10208650
11966406
8990122
8791781
10669647
10785505
11273991
11791158
11983026
12042665
12042666
12042671
12142726
12142729
12172212
12172217
12172220
12194916
12235453
12270762
12360106
12374873
12439220
12439223
12471611
12490595
12544511
12608943
12618594
12618863
12619042
12624279
12668920
12668921
12709726
12724616
12724617
12724621
12724623
12746735
12777962
12781330
12782962
12795791
12798882
12801479
20004378
20087150
20100220
20111011
20132287
20155948
20169140
20201778
20202156
20216355
20225210
20250763
20252989
20324596
20339098
20342856
20349854
20366163
20377823
20404014
20417758
20447473
21175742
21201291
92071833
92146398
95166809
95329110
96152237
96192469
96220220
97163946
97243270
97413707
98143523
98180550
98354687
98392843
99229027
99290598
99301892
99308324
"""
