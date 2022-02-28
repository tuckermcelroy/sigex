# code to get calendar date for Sunday of the first week
first.day <- 1
all.date <- weekly2date(first.day,begin,T)
start.date <- all.date[[1]]
end.date <- all.date[[2]]


# ---- * holiday regressors ----

easter.dates <- read.table("data/easter500.txt")
easter.reg <- gethol(easter.dates,7,0,start.date,end.date)

nyd.dates <- read.table("data/newyear500.txt")
nyd.reg <- gethol(nyd.dates,7,0,start.date,end.date)

mlk.dates <- read.table("data/mlk500.txt")
mlk.reg <- gethol(mlk.dates,7,0,start.date,end.date)

gw.dates <- read.table("data/gw500.txt")
gw.reg <- gethol(gw.dates,7,0,start.date,end.date)

mem.dates <- read.table("data/mem500.txt")
mem.reg <- gethol(mem.dates,7,0,start.date,end.date)

ind.dates <- read.table("data/ind500.txt")
ind.reg <- gethol(ind.dates,7,0,start.date,end.date)

labor.dates <- read.table("data/labor500.txt")
labor.reg <- gethol(labor.dates,7,0,start.date,end.date)

col.dates <- read.table("data/columbus500.txt")
col.reg <- gethol(col.dates,7,0,start.date,end.date)

vet.dates <- read.table("data/vet500.txt")
vet.reg <- gethol(vet.dates,7,0,start.date,end.date)

tg.dates <- read.table("data/thanksgiving500.txt")
tg.reg <- gethol(tg.dates,7,0,start.date,end.date)

xmas.dates <- read.table("data/xmas500.txt")
xmas.reg <- gethol(xmas.dates,7,0,start.date,end.date)

black.dates <- read.table("data/black400.txt")
black.reg <- gethol(black.dates,7,0,start.date,end.date)

## Independence Day, Veteran's Day, and Christmas are purely seasonal
sum(ind.reg^2)
sum(vet.reg^2)
sum(xmas.reg^2)

# ---- * weekly flow regressors ----

easter.reg <- sigex.daily2weekly(easter.reg,first.day,start.date)
easter.reg <- rowSums(easter.reg)/7

nyd.reg <- sigex.daily2weekly(nyd.reg,first.day,start.date)
nyd.reg <- rowSums(nyd.reg)/7

mlk.reg <- sigex.daily2weekly(mlk.reg,first.day,start.date)
mlk.reg <- rowSums(mlk.reg)/7

gw.reg <- sigex.daily2weekly(gw.reg,first.day,start.date)
gw.reg <- rowSums(gw.reg)/7

mem.reg <- sigex.daily2weekly(mem.reg,first.day,start.date)
mem.reg <- rowSums(mem.reg)/7

ind.reg <- sigex.daily2weekly(ind.reg,first.day,start.date)
ind.reg <- rowSums(ind.reg)/7

labor.reg <- sigex.daily2weekly(labor.reg,first.day,start.date)
labor.reg <- rowSums(labor.reg)/7

col.reg <- sigex.daily2weekly(col.reg,first.day,start.date)
col.reg <- rowSums(col.reg)/7

vet.reg <- sigex.daily2weekly(vet.reg,first.day,start.date)
vet.reg <- rowSums(vet.reg)/7

tg.reg <- sigex.daily2weekly(tg.reg,first.day,start.date)
tg.reg <- rowSums(tg.reg)/7

xmas.reg <- sigex.daily2weekly(xmas.reg,first.day,start.date)
xmas.reg <- rowSums(xmas.reg)/7

black.reg <- sigex.daily2weekly(black.reg,first.day,start.date)
black.reg <- rowSums(black.reg)/7
