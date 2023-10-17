########################
#### Script for NDC Data
########################

## wipe
rm(list=ls())

library(devtools)
library(Rcpp)

# suppose directory is set to where sigex is located, e.g.
#setwd("C:\\Users\\neide\\Documents\\GitHub\\sigex")
load_all(".")
root.dir <- getwd()
setwd(paste(root.dir,"/tests/NDC",sep=""))

######################
### Part I: load data

# automatic

#############################################################
### Part II: Metadata Specifications and Exploratory Analysis

start.date <- c(1992,1)
end.date <- c(2020,5)
period <- 12

## create ts object and plot
dataALL.ts <- sigex.load(ndc,start.date,period,
                         c("Shipments","NewOrders"),TRUE)


#############################
## select span and transforms

## all data with log transform
transform <- "none"
aggregate <- FALSE
subseries <- c(1,2)
begin.date <- start(dataALL.ts)
end.date <- end(dataALL.ts)
range <- NULL
data.ts <- sigex.prep(dataALL.ts,transform,aggregate,subseries,range,TRUE)


###############################
### Part III: Model Declaration

N <- dim(data.ts)[2]
T <- dim(data.ts)[1]

###################
## Basic Model: VAR

## preliminary analysis
ar.fit <- ar.yw(diff(ts(ndc[2:T,])))
p.order <- ar.fit$order
par.yw <- aperm(ar.fit$ar,c(2,3,1))
covmat.yw <- getGCD(ar.fit$var.pred,2)
var.out <- var2.par2pre(par.yw)
psi.init <- as.vector(c(covmat.yw[[1]][2,1],log(covmat.yw[[2]]),
                        var.out,colMeans(diff(ts(ndc[2:T,])))))

## model construction
mdl <- NULL
mdl <- sigex.add(mdl,seq(1,N),"varma",c(p.order,0),NULL,"process",c(1,-1))
# regressors:
mdl <- sigex.meaninit(mdl,data.ts,0)


#############################
### Part IV: Model Estimation

## parameter initialization
constraint <- NULL
psi.mle <- psi.init
par.mle <- sigex.psi2par(psi.mle,mdl,data.ts)

## run fitting:
fit.mle <- sigex.mlefit(data.ts,par.mle,constraint,mdl,"bfgs",debug=TRUE)

# fit.mle <-
#   list(list(par = c(1.53390028981385, 14.242031592804, 17.0324362003599,
#                     -0.471755559540316, 0.48972706925505, -0.0176770274069112, -0.392498776720579,
#                     -0.399094699431813, -0.0521993907276512, 0.161306469612418, -0.300580408212516,
#                     0.173188354719214, 0.36650806964924, 0.181795873784704, -0.064448529406187,
#                     0.136113266131282, 0.496532823253715, 0.013062010813688, -0.0818714207759865,
#                     0.0102111906250288, 0.445693306519572, 0.0898406502925605, -0.268826665914594,
#                     0.0903360821499277, 0.000320637689799219, 0.0484855128438093,
#                     -0.144308675343334, -0.0089379765503397, -0.0907715841888581,
#                     0.0516024410959549, -0.0326240049057978, 75.9102025755482, 77.8748318732065
#   ), value = 11290.1700235853, counts = c(`function` = 77L, gradient = 77L
#   ), convergence = 0L, message = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH",
#   hessian = structure(c(41.7670239585277, 0.066904703999171,
#                         0.15543878362223, -0.242634769165306, 0.0905547494767234,
#                         -0.104534365164, 0.35585367186286, 0.225461462832754, 0.170873136084992,
#                         -1.13159489956161, 0.0432885371992597, 0.45453248276317,
#                         0.0148784238263033, 0.209234258363722, 0.0970812834566459,
#                         -0.336737457473646, 0.159565843205201, -0.157904651132412,
#                         0.268451685769833, 0.819444267108338, 0.0105119397630915,
#                         -0.332600166075281, 0.159975797942025, 0.731179852664354,
#                         -0.0208867731998907, 0.640164216747507, 0.111204826680478,
#                         0.368158453056822, -0.0803236162028043, 0.0902532519830856,
#                         0.130194621306146, -0.00321915649692528, 0.00104341779660899,
#                         0.066904703999171, 339.521233399864, -0.371527676179539,
#                         -1.02913145383354, 0.0792078935774043, -0.267630184680456,
#                         0.375539912056411, -0.899432734513539, 0.275610545941163,
#                         -1.36930202643271, -0.270280679615098, 1.80815550265834,
#                         -0.421940967498813, 2.71673616225598, -0.0295360678137513,
#                         -0.395384631701745, 0.0735271896701306, 0.512565293320222,
#                         -0.0185639237315627, 0.976545152298058, -0.00192403604160063,
#                         0.606983121542726, 0.128102101371041, 1.97213239516714, -0.212540726352017,
#                         3.40343535754073, -0.4501489456743, -0.858102112033521, 0.0845791419123998,
#                         1.5278262708307, -0.0261982222582446, -0.00202453520614654,
#                         0.00118507159641013, 0.15543878362223, -0.371527676179539,
#                         339.782259288768, 0.45291631067812, -0.849612661113497, 0.474822172691347,
#                         -0.934725449042162, -0.244000148086343, 0.264178652287228,
#                         2.0259763005015, -0.746514615457272, -0.639737891106051,
#                         0.0860804902913515, 2.65123412646062, -0.692620005793287,
#                         -0.277563458439545, -0.137880761030829, -2.52420045399049,
#                         -0.311473513647798, 0.0863294644659618, -0.284455836663255,
#                         3.48949447470659, -1.00482657217071, -3.79510152015428, 1.18972980089893,
#                         1.46989918903273, 0.222608150579617, -2.65646417574317, 0.556569148102426,
#                         6.28437487648625, -0.616038050793577, 0.0153586370288394,
#                         -0.0114598606160143, -0.242634769165306, -1.02913145383354,
#                         0.45291631067812, 1521.11140641864, 304.674918425007, 1481.69353587946,
#                         -469.443363726896, -908.530464585056, -388.968447623483,
#                         733.283296085574, -104.549504612805, 99.1676070043468, -199.057364625332,
#                         -242.425405758695, 29.6634834739962, 102.453801446245, -163.196803896426,
#                         429.538317348488, -225.676012632903, -1146.31795941023, -273.87741874918,
#                         -364.398458714277, -2.86198746834998, -555.61071962984, 107.161985397397,
#                         556.049475108011, -98.3876316240639, -11.6580818030343, 55.7654211661429,
#                         660.795786188828, -19.4387246210681, -0.00643603925709613,
#                         0.00571230884816032, 0.0905547494767234, 0.0792078935774043,
#                         -0.849612661113497, 304.674918425007, 967.814871273731, 2126.03552290602,
#                         116.355911359278, -1051.24050242011, -211.385357943072, -1549.66284571856,
#                         42.1524450757715, 777.690744143911, -259.277434452088, -539.704658422124,
#                         171.888029854017, -314.488367621379, 140.485415158764, 590.251886933402,
#                         -138.155117838323, -472.203692879702, -218.920862607774,
#                         -804.277283805277, 250.565922669921, 673.800170716277, -44.2975083387864,
#                         312.650605337694, 8.89452849150985, 725.261097613839, 8.00981115389732,
#                         -316.301670409302, 71.5731177933776, -0.00393129084841348,
#                         -8.07176547823474e-05, -0.104534365164, -0.267630184680456,
#                         0.474822172691347, 1481.69353587946, 2126.03552290602, 7841.28836221498,
#                         -990.983458905248, -3355.60184021233, -1127.11655219755,
#                         -1444.54143901385, -246.631818299647, 1962.20444968276, -971.693268184026,
#                         -356.299056875287, 618.537752416159, -1190.34805629781, 196.6624099623,
#                         115.439160708775, -499.176735957008, -2230.55042351916, -1159.70816887057,
#                         -395.460056097363, 312.14922091749, -1350.01532271417, 280.784342521656,
#                         1048.24832851591, -277.312932666973, 2038.87782026868, -48.3266355786327,
#                         472.758185878774, 109.595627236558, -0.0293100583803607,
#                         0.0153941073222086, 0.35585367186286, 0.375539912056411,
#                         -0.934725449042162, -469.443363726896, 116.355911359278,
#                         -990.983458905248, 1027.67309363117, 357.226211690431, 261.466486790596,
#                         -1515.48157964498, 118.551518880849, 129.309752537665, -36.2395680895133,
#                         660.132768643962, -76.6531361477973, -407.176901717321, 26.8398976004391,
#                         -222.05643426787, 251.441676937247, 826.899570711248, 140.874379212619,
#                         -205.743635433464, 91.3749640858441, 876.467113812396, -223.313842070638,
#                         -417.480156556849, 123.08016539464, 20.1361904146324, 12.4543548736256,
#                         -397.443636757089, 45.5322228845034, 0.0153822838910855,
#                         -0.0141835698741488, 0.225461462832754, -0.899432734513539,
#                         -0.244000148086343, -908.530464585056, -1051.24050242011,
#                         -3355.60184021233, 357.226211690431, 1919.61352402359, 388.968297556858,
#                         998.813817659538, 51.1372172695701, -835.422732052393, 372.666550902068,
#                         387.967514598131, -329.222486925573, 473.23987882919, -155.118916609354,
#                         -14.8546168929897, 178.281291482563, 864.867993186635, 467.592008590145,
#                         443.9849026312, -179.372042566683, 178.722293185274, -86.133250988496,
#                         -591.482325489778, 121.964692198162, -1053.95097398286, 46.0951118839148,
#                         128.3021956624, -75.252238275425, 0.0311129042529501, -0.0143882061820477,
#                         0.170873136084992, 0.275610545941163, 0.264178652287228,
#                         -388.968447623483, -211.385357943072, -1127.11655219755,
#                         261.466486790596, 388.968297556858, 558.927455585945, -1004.63566104736,
#                         -43.9059649579576, 46.083449433354, 161.17116138048, 619.113128777826,
#                         -209.828536753776, -229.322338782367, 252.389599154412, -554.907919195102,
#                         118.226556878653, 905.744919691642, 190.72402801612, 93.2246912270784,
#                         72.1012866051751, 559.322483695723, 5.47425634067622, -122.710282084881,
#                         13.7263814394828, -93.4973450057441, 25.2658173849341, -334.966548052762,
#                         -32.6758931805671, -0.00223076312977355, 0.00401109900849406,
#                         -1.13159489956161, -1.36930202643271, 2.0259763005015, 733.283296085574,
#                         -1549.66284571856, -1444.54143901385, -1515.48157964498,
#                         998.813817659538, -1004.63566104736, 8264.73506594994, -266.337501670932,
#                         -2578.95066579295, 627.302244993189, -1401.41121164561, -220.336103211594,
#                         1763.28451152585, -985.618238246389, -201.264148017799, -152.45954546117,
#                         -1967.81691374781, -43.8793476860155, 81.3940205262043, -434.335639056371,
#                         -3184.63530356894, 331.288580127875, 217.441141558083, -51.3945333295851,
#                         -1260.67726341716, 112.98723939035, 1663.635055138, -81.6432211649953,
#                         -0.0409613676310983, 0.0221969003177946, 0.0432885371992597,
#                         -0.270280679615098, -0.746514615457272, -104.549504612805,
#                         42.1524450757715, -246.631818299647, 118.551518880849, 51.1372172695701,
#                         -43.9059649579576, -266.337501670932, 413.678089898895, -83.9127778817783,
#                         -118.835816465435, -434.903159202804, 189.722439245088, 131.333571971481,
#                         118.223914341797, -30.6096571875969, 123.240111861378, 187.283284503792,
#                         98.9633249446342, -252.887227816245, -2.34626304518315, 538.442408696937,
#                         -214.577466522314, -162.694845130318, 40.5773798775044, 257.730505381915,
#                         -158.823390847829, -150.735920215084, 82.9631039778178, 0.00788759280112572,
#                         -0.00945647116168402, 0.45453248276317, 1.80815550265834,
#                         -0.639737891106051, 99.1676070043468, 777.690744143911, 1962.20444968276,
#                         129.309752537665, -835.422732052393, 46.083449433354, -2578.95066579295,
#                         -83.9127778817783, 1797.11564123863, -553.691543700552, 762.385110192554,
#                         62.2874049440725, -863.46199577747, 264.018646703335, 416.243311065045,
#                         -167.400036389154, 0.117322997539304, -284.157913029048,
#                         200.204461634712, 126.361934235319, 287.526567717578, 18.3818485766096,
#                         550.206357274874, -100.713092251681, 240.247006786376, 34.0095398314588,
#                         -308.28111766823, 20.8692638352659, 0.0428335624746978, -0.0166944573720684,
#                         0.0148784238263033, -0.421940967498813, 0.0860804902913515,
#                         -199.057364625332, -259.277434452088, -971.693268184026,
#                         -36.2395680895133, 372.666550902068, 161.17116138048, 627.302244993189,
#                         -118.835816465435, -553.691543700552, 513.188180775614, -964.296846177604,
#                         3.25623432217981, 561.914238005556, -55.0041013411828, -44.5959558419418,
#                         137.943959089171, 160.19710483306, 202.767897008016, -184.696184533095,
#                         30.9595980070299, 461.676451777748, 53.718449180451, -312.087291604257,
#                         100.22359992945, 179.717430455639, -10.3686747934262, -143.273906815011,
#                         -48.0986705042596, -0.0159111550601665, 0.00911177266971208,
#                         0.209234258363722, 2.71673616225598, 2.65123412646062, -242.425405758695,
#                         -539.704658422124, -356.299056875287, 660.132768643962, 387.967514598131,
#                         619.113128777826, -1401.41121164561, -434.903159202804, 762.385110192554,
#                         -964.296846177604, 8059.22949757587, -923.358153158915, -3377.3218819988,
#                         255.846489835676, -1646.29743812839, -382.514242119214, 707.853012045234,
#                         -282.154961723791, 267.521089426737, -179.614846729237, -2838.49925881441,
#                         201.024522993976, 346.952724157745, -189.194283848337, -2531.73956571118,
#                         217.245022213319, 652.21471913901, -70.7725821484928, -0.0176717094291234,
#                         0.0125714905152563, 0.0970812834566459, -0.0295360678137513,
#                         -0.692620005793287, 29.6634834739962, 171.888029854017, 618.537752416159,
#                         -76.6531361477973, -329.222486925573, -209.828536753776,
#                         -220.336103211594, 189.722439245088, 62.2874049440725, 3.25623432217981,
#                         -923.358153158915, 609.444160545536, 152.720348523872, -24.7270863837912,
#                         -267.424583398679, 113.586404950183, -7.4389745350345, 14.1822583827889,
#                         57.5048884456919, 25.1737978942401, 373.983394638344, -175.071754028977,
#                         -21.5491104427201, -9.27462451727479, 740.838150704803, -237.949017900974,
#                         -142.198812682182, 58.8102657275158, 0.00553473000763915,
#                         -0.00801765054347925, -0.336737457473646, -0.395384631701745,
#                         -0.277563458439545, 102.453801446245, -314.488367621379,
#                         -1190.34805629781, -407.176901717321, 473.23987882919, -229.322338782367,
#                         1763.28451152585, 131.333571971481, -863.46199577747, 561.914238005556,
#                         -3377.3218819988, 152.720348523872, 2569.09904555869, -373.622840925236,
#                         1050.11504638242, 120.920552490134, -249.458588314155, 255.106350778078,
#                         554.567400286032, -84.0422580949962, 866.388552367425, -43.055113110313,
#                         -307.409489323618, 126.508275570814, 384.027151540067, 0.472934516437817,
#                         391.168477563042, -42.0041556026263, 0.0360123522114009,
#                         -0.0140905740408925, 0.159565843205201, 0.0735271896701306,
#                         -0.137880761030829, -163.196803896426, 140.485415158764,
#                         196.6624099623, 26.8398976004391, -155.118916609354, 252.389599154412,
#                         -985.618238246389, 118.223914341797, 264.018646703335, -55.0041013411828,
#                         255.846489835676, -24.7270863837912, -373.622840925236, 373.400216176378,
#                         -368.258179150871, 78.4060389378283, 354.990839241509, 30.7155019072525,
#                         -127.697747302591, 106.651621081255, 356.476866727462, -18.8424619409489,
#                         -192.514572518121, 1.61996831593569, 354.974976289668, -78.5579686635174,
#                         -344.543726441771, 46.813940116408, -0.0105137587524951,
#                         0.00469071892439388, -0.157904651132412, 0.512565293320222,
#                         -2.52420045399049, 429.538317348488, 590.251886933402, 115.439160708775,
#                         -222.05643426787, -14.8546168929897, -554.907919195102, -201.264148017799,
#                         -30.6096571875969, 416.243311065045, -44.5959558419418, -1646.29743812839,
#                         -267.424583398679, 1050.11504638242, -368.258179150871, 8164.67074037064,
#                         -731.032630028494, -3167.31150269334, 145.489762417128, -1661.69052226905,
#                         63.7176785858173, 782.196780619415, -83.5097807794227, 209.007473131351,
#                         -13.801500244881, -2724.27204936321, 240.593693888513, 134.522519601887,
#                         -6.00613657297799, -0.00901854946278036, 0.0075724528869614,
#                         0.268451685769833, -0.0185639237315627, -0.311473513647798,
#                         -225.676012632903, -138.155117838323, -499.176735957008,
#                         251.441676937247, 178.281291482563, 118.226556878653, -152.45954546117,
#                         123.240111861378, -167.400036389154, 137.943959089171, -382.514242119214,
#                         113.586404950183, 120.920552490134, 78.4060389378283, -731.032630028494,
#                         539.174634013762, 494.771379635495, -14.9937036439951, 146.227451295999,
#                         -51.4740263497515, 444.518876065558, -83.6548088045674, -50.4753725181217,
#                         10.0083402685414, 431.685520197789, -179.738524366257, -173.436564182339,
#                         15.8017978719727, 0.0141694727062713, -0.0119644028018229,
#                         0.819444267108338, 0.976545152298058, 0.0863294644659618,
#                         -1146.31795941023, -472.203692879702, -2230.55042351916,
#                         826.899570711248, 864.867993186635, 905.744919691642, -1967.81691374781,
#                         187.283284503792, 0.117322997539304, 160.19710483306, 707.853012045234,
#                         -7.4389745350345, -249.458588314155, 354.990839241509, -3167.31150269334,
#                         494.771379635495, 3692.95046539264, 291.892858967913, 1318.96596985825,
#                         -10.0077586466796, 1346.73519914941, -162.416606144689, -188.734428775206,
#                         147.328574485073, 994.703248125006, -90.196283053956, -833.888867418864,
#                         16.0553768182581, 0.0557236035092501, -0.0258710315392818,
#                         0.0105119397630915, -0.00192403604160063, -0.284455836663255,
#                         -273.87741874918, -218.920862607774, -1159.70816887057, 140.874379212619,
#                         467.592008590145, 190.72402801612, -43.8793476860155, 98.9633249446342,
#                         -284.157913029048, 202.767897008016, -282.154961723791, 14.1822583827889,
#                         255.106350778078, 30.7155019072525, 145.489762417128, -14.9937036439951,
#                         291.892858967913, 319.676815706771, -167.253455401806, 56.932114148367,
#                         531.801304532564, -98.3022068794526, -390.307530778955, 91.3516287255334,
#                         -106.274090285297, -17.4172773768078, -239.638846323942,
#                         -1.06291872725706, -0.00343788997270167, 0.000766704033594579,
#                         -0.332600166075281, 0.606983121542726, 3.48949447470659,
#                         -364.398458714277, -804.277283805277, -395.460056097363,
#                         -205.743635433464, 443.9849026312, 93.2246912270784, 81.3940205262043,
#                         -252.887227816245, 200.204461634712, -184.696184533095, 267.521089426737,
#                         57.5048884456919, 554.567400286032, -127.697747302591, -1661.69052226905,
#                         146.227451295999, 1318.96596985825, -167.253455401806, 7932.16215652137,
#                         -662.57121397939, -3169.82528738663, 240.086028497899, -1755.46985747133,
#                         120.518438052386, 533.805900886364, -52.1631432093272, 142.744977893017,
#                         -10.6926058833778, -0.00463137439510319, 0.00498448571306653,
#                         0.159975797942025, 0.128102101371041, -1.00482657217071,
#                         -2.86198746834998, 250.565922669921, 312.14922091749, 91.3749640858441,
#                         -179.372042566683, 72.1012866051751, -434.335639056371, -2.34626304518315,
#                         126.361934235319, 30.9595980070299, -179.614846729237, 25.1737978942401,
#                         -84.0422580949962, 106.651621081255, 63.7176785858173, -51.4740263497515,
#                         -10.0077586466796, 56.932114148367, -662.57121397939, 445.656801730365,
#                         503.146110077068, -174.607134795224, 246.366214014415, -88.906369455799,
#                         181.08534231942, 33.4074375132332, -120.605001029617, 1.81423138201353,
#                         0.00896443452802487, -0.00692375579092186, 0.731179852664354,
#                         1.97213239516714, -3.79510152015428, -555.61071962984, 673.800170716277,
#                         -1350.01532271417, 876.467113812396, 178.722293185274, 559.322483695723,
#                         -3184.63530356894, 538.442408696937, 287.526567717578, 461.676451777748,
#                         -2838.49925881441, 373.983394638344, 866.388552367425, 356.476866727462,
#                         782.196780619415, 444.518876065558, 1346.73519914941, 531.801304532564,
#                         -3169.82528738663, 503.146110077068, 5870.19973318093, -563.329531360068,
#                         305.898342503497, 157.53037041577, 1488.95452139186, -161.105926508753,
#                         -751.343049159914, 1.7946690604731, 0.0917821125767659, -0.0458953763882164,
#                         -0.0208867731998907, -0.212540726352017, 1.18972980089893,
#                         107.161985397397, -44.2975083387864, 280.784342521656, -223.313842070638,
#                         -86.133250988496, 5.47425634067622, 331.288580127875, -214.577466522314,
#                         18.3818485766096, 53.718449180451, 201.024522993976, -175.071754028977,
#                         -43.055113110313, -18.8424619409489, -83.5097807794227, -83.6548088045674,
#                         -162.416606144689, -98.3022068794526, 240.086028497899, -174.607134795224,
#                         -563.329531360068, 307.772421365371, -22.305624497676, 24.1666712099686,
#                         -92.255967501842, 100.016673968639, 90.8158799575176, -18.4266959877277,
#                         -0.0255035956797656, 0.0178645223058993, 0.640164216747507,
#                         3.40343535754073, 1.46989918903273, 556.049475108011, 312.650605337694,
#                         1048.24832851591, -417.480156556849, -591.482325489778, -122.710282084881,
#                         217.441141558083, -162.694845130318, 550.206357274874, -312.087291604257,
#                         346.952724157745, -21.5491104427201, -307.409489323618, -192.514572518121,
#                         209.007473131351, -50.4753725181217, -188.734428775206, -390.307530778955,
#                         -1755.46985747133, 246.366214014415, 305.898342503497, -22.305624497676,
#                         8488.77262069436, -780.022078288312, -3325.05905680591, 310.733021706255,
#                         -1535.91825664989, 140.255396217981, -0.0230488694796804,
#                         0.0153536348079797, 0.111204826680478, -0.4501489456743,
#                         0.222608150579617, -98.3876316240639, 8.89452849150985, -277.312932666973,
#                         123.08016539464, 121.964692198162, 13.7263814394828, -51.3945333295851,
#                         40.5773798775044, -100.713092251681, 100.22359992945, -189.194283848337,
#                         -9.27462451727479, 126.508275570814, 1.61996831593569, -13.801500244881,
#                         10.0083402685414, 147.328574485073, 91.3516287255334, 120.518438052386,
#                         -88.906369455799, 157.53037041577, 24.1666712099686, -780.022078288312,
#                         481.500322621287, 316.398473842128, -177.801809968514, 85.8534549479373,
#                         -90.1356479516835, 0.0114878275780939, -0.00958721102506388,
#                         0.368158453056822, -0.858102112033521, -2.65646417574317,
#                         -11.6580818030343, 725.261097613839, 2038.87782026868, 20.1361904146324,
#                         -1053.95097398286, -93.4973450057441, -1260.67726341716,
#                         257.730505381915, 240.247006786376, 179.717430455639, -2531.73956571118,
#                         740.838150704803, 384.027151540067, 354.974976289668, -2724.27204936321,
#                         431.685520197789, 994.703248125006, -106.274090285297, 533.805900886364,
#                         181.08534231942, 1488.95452139186, -92.255967501842, -3325.05905680591,
#                         316.398473842128, 5622.11033366111, -480.45136327346, 117.566962217097,
#                         14.7921105053683, 0.0687823558109812, -0.0357736098521855,
#                         -0.0803236162028043, 0.0845791419123998, 0.556569148102426,
#                         55.7654211661429, 8.00981115389732, -48.3266355786327, 12.4543548736256,
#                         46.0951118839148, 25.2658173849341, 112.98723939035, -158.823390847829,
#                         34.0095398314588, -10.3686747934262, 217.245022213319, -237.949017900974,
#                         0.472934516437817, -78.5579686635174, 240.593693888513, -179.738524366257,
#                         -90.196283053956, -17.4172773768078, -52.1631432093272, 33.4074375132332,
#                         -161.105926508753, 100.016673968639, 310.733021706255, -177.801809968514,
#                         -480.45136327346, 320.615854434436, 29.3005427920434, 4.91519267598051,
#                         -0.0201162038138136, 0.0153972905536648, 0.0902532519830856,
#                         1.5278262708307, 6.28437487648625, 660.795786188828, -316.301670409302,
#                         472.758185878774, -397.443636757089, 128.3021956624, -334.966548052762,
#                         1663.635055138, -150.735920215084, -308.28111766823, -143.273906815011,
#                         652.21471913901, -142.198812682182, 391.168477563042, -344.543726441771,
#                         134.522519601887, -173.436564182339, -833.888867418864, -239.638846323942,
#                         142.744977893017, -120.605001029617, -751.343049159914, 90.8158799575176,
#                         -1535.91825664989, 85.8534549479373, 117.566962217097, 29.3005427920434,
#                         8706.21140211369, -685.885469465575, -0.0273678324447246,
#                         0.017943875718629, 0.130194621306146, -0.0261982222582446,
#                         -0.616038050793577, -19.4387246210681, 71.5731177933776,
#                         109.595627236558, 45.5322228845034, -75.252238275425, -32.6758931805671,
#                         -81.6432211649953, 82.9631039778178, 20.8692638352659, -48.0986705042596,
#                         -70.7725821484928, 58.8102657275158, -42.0041556026263, 46.813940116408,
#                         -6.00613657297799, 15.8017978719727, 16.0553768182581, -1.06291872725706,
#                         -10.6926058833778, 1.81423138201353, 1.7946690604731, -18.4266959877277,
#                         140.255396217981, -90.1356479516835, 14.7921105053683, 4.91519267598051,
#                         -685.885469465575, 462.220623830945, 0.00882755557540804,
#                         -0.00953855305851903, -0.00321915649692528, -0.00202453520614654,
#                         0.0153586370288394, -0.00643603925709613, -0.00393129084841348,
#                         -0.0293100583803607, 0.0153822838910855, 0.0311129042529501,
#                         -0.00223076312977355, -0.0409613676310983, 0.00788759280112572,
#                         0.0428335624746978, -0.0159111550601665, -0.0176717094291234,
#                         0.00553473000763915, 0.0360123522114009, -0.0105137587524951,
#                         -0.00901854946278036, 0.0141694727062713, 0.0557236035092501,
#                         -0.00343788997270167, -0.00463137439510319, 0.00896443452802487,
#                         0.0917821125767659, -0.0255035956797656, -0.0230488694796804,
#                         0.0114878275780939, 0.0687823558109812, -0.0201162038138136,
#                         -0.0273678324447246, 0.00882755557540804, 0.00151248968904838,
#                         -0.00081672624219209, 0.00104341779660899, 0.00118507159641013,
#                         -0.0114598606160143, 0.00571230884816032, -8.07176547823474e-05,
#                         0.0153941073222086, -0.0141835698741488, -0.0143882061820477,
#                         0.00401109900849406, 0.0221969003177946, -0.00945647116168402,
#                         -0.0166944573720684, 0.00911177266971208, 0.0125714905152563,
#                         -0.00801765054347925, -0.0140905740408925, 0.00469071892439388,
#                         0.0075724528869614, -0.0119644028018229, -0.0258710315392818,
#                         0.000766704033594579, 0.00498448571306653, -0.00692375579092186,
#                         -0.0458953763882164, 0.0178645223058993, 0.0153536348079797,
#                         -0.00958721102506388, -0.0357736098521855, 0.0153972905536648,
#                         0.017943875718629, -0.00953855305851903, -0.00081672624219209,
#                         0.000506588548887521), .Dim = c(33L, 33L))), list(list(structure(c(1,
#                                                                                            1.53390028981385, 0, 1), .Dim = c(2L, 2L))), list(c(14.242031592804,
#                                                                                                                                                17.0324362003599)), list(structure(c(-0.515491359942004, 0.402577194696,
#                                                                                                                                                                                     0.0922318695603973, -0.625048092872619, -0.298176005683577, 0.39498073180832,
#                                                                                                                                                                                     0.095863361014835, -0.52322291642085, 0.187895825290406, 0.901710614761734,
#                                                                                                                                                                                     0.0944632126558633, -0.334184976053577, 0.151677316698701, 0.925120123908555,
#                                                                                                                                                                                     0.0393623104141371, -0.302800180493954, 0.0412911410046517, 0.527745200813104,
#                                                                                                                                                                                     0.0669552957585035, -0.25151116668032, 0.081131610454198, -0.0575125801237395,
#                                                                                                                                                                                     0.0532375307220312, -0.0877465918304333, -0.0154166849609989,
#                                                                                                                                                                                     -0.113054154478938, 0.041221695164254, -0.0185242254002529), .Dim = c(2L,
#                                                                                                                                                                                                                                                           2L, 7L))), c(75.9102025755482, 77.8748318732065)))

## MLE fitting results
#  divergence:    11290.17
# psi.mle <- c(1.53390034116153, 14.2420315924048, 17.0324361912539, -0.471755451431251,
# 0.489727692311396, -0.017677234835485, -0.392499692057969, -0.39909419297684,
# -0.0522002118651599, 0.161306450376942, -0.300580812314135, 0.173188276541793,
# 0.366507618868933, 0.181796113107285, -0.0644480542794434, 0.136113149623089,
# 0.496532863295004, 0.013062477268366, -0.0818707058555948, 0.0102111694634749,
# 0.445693008411153, 0.0898405276107301, -0.268826276410836, 0.0903359864828166,
# 0.000320560369556658, 0.0484856870030809, -0.144308229871243,
# -0.0089375652660538, -0.0907709494584039, 0.0516025826426896,
# -0.0326238274507487, 75.9102025269533, 77.8748319106466)

## manage output
psi.mle <- sigex.eta2psi(fit.mle[[1]]$par,constraint)
hess <- fit.mle[[1]]$hessian
par.mle <- fit.mle[[2]]

## residual analysis
resid.mle <- sigex.resid(psi.mle,mdl,data.ts)[[1]]
resid.mle <- sigex.load(t(resid.mle),start(data.ts),frequency(data.ts),
                        colnames(data.ts),TRUE)
resid.acf <- acf(resid.mle,lag.max=4*period,plot=TRUE)$acf

## examine condition numbers
log(sigex.conditions(data.ts,psi.mle,mdl))

## model checking
sigex.portmanteau(resid.mle,4*period,length(psi.mle))
sigex.gausscheck(resid.mle)

## check on standard errors and get t statistics
print(eigen(hess)$values)
tstats <- sigex.tstats(mdl,psi.mle,hess,constraint)
print(tstats)

## bundle
analysis.mle <- sigex.bundle(data.ts,transform,mdl,psi.mle)


##########################################
### Part V: Casting

## load up the fitted model for casting
data.ts <- analysis.mle[[1]]
mdl <- analysis.mle[[3]]
psi <- analysis.mle[[4]]
param <- sigex.psi2par(psi,mdl,data.ts)

## Generate aftcasts and forecasts with uncertainty
window <- 50
data.casts <- sigex.midcast(psi,mdl,data.ts,window)
extract.casts <- sigex.castextract(data.ts,data.casts,mdl,window,param)

## display
castcol <- "black"
fade <- 60
dataPad.ts <- rbind(matrix(NA,nrow=window,ncol=N),data.ts,matrix(NA,nrow=window,ncol=N))
#pdf(file="NdcCasts.pdf",height=8,width=10)
par(mfrow=c(2,1))
for(i in 1:N)
{
  plot(ts(dataPad.ts[,i],start=start.date,frequency=period),
       xlab="Year",ylab="",lwd=1,col=1)
  sigex.graph(extract.casts,NULL,start.date,period,i,0,castcol,fade)
}
dev.off()

