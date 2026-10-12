|Population                        |CLT assumptions                |Predicted n |Observed n |Skew @ n=30 |SD/SE @ n=30 |Tail rate @ n=30 |Normal @ n=30? |
|:---------------------------------|:------------------------------|:-----------|:----------|:-----------|:------------|:----------------|:--------------|
|Standard Normal                   |None (i.i.d., finite variance) |1           |2          |-0.01       |1.01         |5.2%             |Yes            |
|Six-Sided Die                     |None (i.i.d., finite variance) |5           |10         |0.00        |1.00         |4.8%             |Yes            |
|Exponential (rate = 1)            |None (i.i.d., finite variance) |100         |100        |0.40        |0.99         |4.8%             |No             |
|Binomial (5, 0.10)                |None (i.i.d., finite variance) |36          |50         |0.21        |1.00         |4.4%             |No             |
|Cauchy (0, 1)                     |Breaks: finite variance        |n/a         |> 2000     |undefined   |undefined    |undefined        |No             |
|Dependent Machine Failure         |Breaks: independent            |n/a         |1000       |0.75        |1.01         |4.4%             |No             |
|Non-identical Bernoulli (Fatigue) |Breaks: identical              |15          |100        |-0.14       |0.99         |5.9%             |Yes            |
|Age at Death (Mixture)            |None (i.i.d., finite variance) |19          |50         |-0.20       |1.01         |5.1%             |No             |
