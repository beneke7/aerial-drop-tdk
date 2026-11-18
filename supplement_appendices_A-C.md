# Függelékek A–C, bővebb munkaváltozat

A TDK-dolgozat (*A légi tűzoltás fizikai modellezése és hatékonyságának vizsgálata*, ELTE TTK, 2026) A–C. függelékének rövidítés előtti változata, 2026. november 18-i állapot. A beadott dolgozatban ezek tömörítve szerepelnek; az itt szereplő kiegészítő táblázatok és megjegyzések a számítások részleteit dokumentálják. Az ábra- és táblázatszámok a végleges dolgozatétól kissé eltérhetnek. A hivatkozások kulcsai a dolgozat irodalomjegyzékére utalnak.

# A. függelék: Kiegészítések a 2. fejezethez

**A folytonos hatásfok.** Szükségletértékenként minden keresztirányú cellaoszlop (2.2. szakasz) a tűz megállításának valószínűségével számít be, a lángok véletlen elérési távolságára átlagolva. $\eta_B$ és $s_B$ a hasznosuló hányad tíz szükségletértékre vett átlaga és szórása.

| Elem | Érték | Forrás / megjegyzés |
|---|---|---|
| A lángok elérési távolsága, nedvesítetlen éghető anyag | logisztikus eloszlás; az esetek felében 1.5, körülbelül 2%-ában 2.5 lánghosszúságnál nagyobb | feltevés; összevetésül az áttörési modell [Wilson1988] szavannai lánghosszúságokkal [Brou2022, 7. o., 18–19. egyenlet] 3–16 MW/m-re 31–51%-ot, illetve 1–6.5%-ot ad |
| Lánghosszúság, részben átnedvesedett éghető anyag | a nedvesítetlen $\beta^{0.46}$-szorosa | $L_f\propto I^{0.46}$ [Byram1959], változatlan $H_c w_c$ (saját származtatás); $\beta$ a terjedési sebességek aránya a Rothermel-modell nedvességi csillapítási tényezőjéből és előmelegítési hőjéből [Rothermel1972, 7. és 12. o.] |
| Nedvességtartalom | nedvesítetlen: a kioltási nedvességtartalom 30%-a; a szükségletnek megfelelő vizet kapott: a kioltási nedvességtartalom | feltevés |
| Szükséglet | a vonal mentén állandó; tíz érték $C^*/1.5$ és $1.5\,C^*$ között, $\ln C$ szerint egyenlő lépések felezőpontjaiban | |

**A geometriai változat.** A 2.4. szakasz harmadik összetevője nélkül a részben átnedvesedett éghető anyagban a tűz úgy terjed és ég, mint a nedvesítetlenben. Oszloponként csak a szükségletet elérő leghosszabb összefüggő szakasz számít be, a szükségletnek pedig 40 értékét vettem.

| Hatásfok-változás vagy érték | Geometriai | Teljes |
|---|---|---|
| Zajmentes egyenletes sáv, sűrű bozót, $1.5\,W C^*$ terhelés, 0.5 / 1 m-es simítás | −13% / −35% | 0 / −14% |
| Az 5. fejezet tartálya, sűrű bozót, simítás 0.5 → 1 m | +31% | +11% |
| Egyenletes / legjobb normál eloszlású sáv a vizsgált terhelések közül a legkedvezőbbnél (2.5. szakasz) | 0.62 / 0.45 | 0.66 / 0.66 |

Zöld bozótban a geometriai változat szerint 0.5 m-es simítással a négy sugár, 0.75 m-es simítással már a tartály kap nagyobb hatásfokot.

*A.1. táblázat.* A kilenc mért szóráskép hatásfoka szintközi interpolációval. Éghető anyagonként a kezelt vonal hossza ($L$, m), a diszkrét ($\eta$) és a folytonos hatásfok ($\eta_B \pm s_B$).

| Kibocsátás | Kibocsátott víz [t] | Legmagasabb feliratozott szint [GPC] | 1. fű | 2. zöld bozót | 3. sűrű bozót |
|---|---|---|---|---|---|
| Black Hawk, Bambi 660, 42 kt* | 2.00 | 10 | 89; 0.06; 0.06 ± 0.00 | 60; 0.17; 0.19 ± 0.01 | 42; 0.18; 0.21 ± 0.04 |
| Black Hawk, Bambi 660, 50 kt* | 2.00 | 4 | 109; 0.08; 0.08 ± 0.00 | 58; 0.17; 0.23 ± 0.07 | n. h. |
| Black Hawk, Bambi 660, 81 kt* | 2.00 | 3 | 102; 0.07; 0.07 ± 0.01 | 0; 0; 0.04 ± 0.04 | n. h. |
| Bell 212, 40 kt | 1.40 | 6.1 | 54; 0.05; 0.05 ± 0.00 | 40; 0.16; 0.17 ± 0.02 | 1; 0.01; 0.06 ± 0.07 |
| S-2T, 118 kt | 4.54 | 3 | 226; 0.07; 0.07 ± 0.01 | 30; 0.04; 0.05 ± 0.05 | n. h. |
| CL-415 | 6.13 | 6.1 | 150; 0.03; 0.03 ± 0.00 | 0; 0; 0.01 ± 0.01 | 0; 0; 0 ± 0 |
| CL-215 | 5.45 | 6.1 | 132; 0.03; 0.03 ± 0.00 | 44; 0.05; 0.06 ± 0.02 | 6; 0.01; 0.01 ± 0.02 |
| AT-802F | 3.03 | 6.1 | 84; 0.04; 0.04 ± 0.00 | 0; 0; 0.03 ± 0.03 | 0; 0; 0 ± 0 |
| DHC-6 | 2.09 | 6.1 | 78; 0.05; 0.05 ± 0.01 | 0; 0; 0.03 ± 0.04 | 0; 0; 0 ± 0 |

n. h.: a szükséges fedettség a legmagasabb feliratozott szint fölött van. Ha e szint a szükséges fedettség 1.5-szeresénél alacsonyabb, $\eta_B$ alsó becslés.
*: feltételezett 2.0 t víz (80%-ig töltött tartály). A 42 kt-s szórásképre sűrű bozótban kapott értékek csak tájékoztató jellegűek, mert 6 GPC-s zárt szintvonalai 239 m²-ből 34 m²-en feliratlanok.

*A.2. táblázat.* A szakirodalomban közölt vízigények l/m²-ben (vízre 1 l/m² = 1 kg/m² = 1 mm).

| Forrás | Fű | Zöld bozót | Sűrű bozót |
|--------|----|----|----|
| Kibocsátási útmutató [Solarz2001, 1. táblázat], átlagos intenzitás, retardáns | 0.41 (1 GPC) | 1.22 (3 GPC) | 2.44 (6 GPC) |
| Francia kézikönyv, Restás [Restas2023, 6. o.] nyomán | 0.4 | 1.6 | 2.4 |
| Kanadai útmutató [McFayden2023, 8-2. táblázat], intenzitás szerint | 0.5 (300–900 kW/m) | 0.5–0.75 (900–1300 kW/m) | 1.5–2.5 (3.2–7.7 MW/m) |
| Áttekintő tanulmány [Wu2026, 14. o.], intenzitás szerint | 0.2 (1.7 MW/m alatt) | 0.2–1.2 (saját becslés) | legalább 1.2 (3 MW/m fölött) |
| Tartomány | 0.2–0.5 | 0.5–1.6 | 1.2–2.5 |
| Legnagyobb / legkisebb | 2.5 | 3.2 | 2.1 |

Az intenzitás szerinti értékeket én rendeltem az éghető anyagokhoz (sűrű bozót: 3–8 MW/m). A zöld bozót tartományából kihagytam a 0.2 l/m²-t, mert vele a legnagyobb és legkisebb érték aránya 8 lenne.

*A.3. táblázat.* Hatásfokok zöld bozótban a simítási hossz szerint. A: teli külső függesztésű tartály, B: az 5. fejezet négy sebességillesztett sugara. Kibocsátásonként négy szóráskép (szembe- és oldalszél, két cseppcsomag-felépítés) átlaga, párolgás nélkül.

| Simítási hossz [m] | $\eta$, A | $\eta$, B | $\eta_B \pm s_B$, A | $\eta_B \pm s_B$, B |
|---|---|---|---|---|
| 0 | 0.017 | 0.042 | 0.094 ± 0.025 | 0.089 ± 0.006 |
| 0.25 | 0.035 | 0.044 | 0.125 ± 0.031 | 0.127 ± 0.008 |
| 0.5 | 0.078 | 0.050 | 0.173 ± 0.031 | 0.158 ± 0.004 |
| 0.75 | 0.113 | 0.125 | 0.195 ± 0.027 | 0.165 ± 0.003 |
| 1 | 0.142 | 0.155 | 0.200 ± 0.025 | 0.167 ± 0.003 |
| 1.5 | 0.171 | 0.159 | 0.194 ± 0.025 | 0.167 ± 0.003 |

**A cellaméret hatása.** Két idealizált sávot, három mért szórásképet és három számított kibocsátást értékeltem ki zöld és sűrű bozótra, változatlan $C^*$-gal és $W$-vel, 0.5, 0.25 és 0.125 m-es cellákon, a két sávot, a három mért szórásképet és a számított kibocsátások szembeszeles változatát 0.1 m-es cellákon is. A mért szórásképek sorrendje egyik hatásfok szerint sem változott, a számítottaké 0.5 m-es simításnál szintén nem.

| Hatásfok, cellaméret | Változás |
|---|---|
| Folytonos, legfeljebb 0.25 m | a legfinomabb rendelkezésre álló rácshoz képest többnyire 2%-on belül, a legnagyobb eltérés 2.85% |
| Folytonos, 0.5 m | többnyire 2%-on belül (a tartály oldalszeles szórásképe, sűrű bozót: 0.058 helyett 0.057); szembeszélben körülbelül 6%-kal kisebb (a tartály, sűrű bozót: 0.284 helyett 0.268) |
| Diszkrét | nem konvergál: zajmentes sávon fél cellányi eltolásra 0 és a teljes érték között ugrik; szórásképeken legfeljebb 10% (Bell 212, sűrű bozót: 67%) |
| A 2.3. ábra három kibocsátásának egy-egy szembeszeles, cellánként egy cseppcsomaggal számított szórásképe sűrű bozótra, 0.25 és 0.125 m-es cellákon | legalább 0.5 m-es simításnál a két rács eredménye a diszkrét hatásfoknál legfeljebb 1.6%-kal, a folytonosnál legfeljebb 2.2%-kal tér el |

**A legkevesebb vizet igénylő normál eloszlású sáv.** A Gauss-alakú keresztmetszet $W$ szélességben eléri a $C^*$ fedettséget, ha

$$
\lambda \ge \sqrt{2\pi}\,C^*\sigma\exp\!\big(W^2/8\sigma^2\big).
$$

A jobb oldal $\sigma = W/2$-nél minimális:

$$
\lambda = \tfrac12\sqrt{2\pi e}\,W C^* = 2.07\,W C^*,\qquad \eta = W C^*/\lambda = 0.48;
$$

ekkor a csúcs $\sqrt{e}\,C^* = 1.65\,C^*$, a sávba a víz $\mathrm{erf}(1/\sqrt2)$ = 68%-a esik.
Nagyobb terhelésnél minden megengedett $\sigma$-ra $\eta = W C^*/\lambda$; a folytonos hatásfok ennek és az oszlopok közös megállítási valószínűségének szorzata.
A 2.5b ábra számítása 0.5 m-es cellákon, 40 szükségletértékkel készült (10 értéknél az egyenletes sáv görbéje $\lambda = W C^*$ közelében 0.07-os lépcsőkkel ugrik); az ábrán a 3.38 m széles egyenletes sáv hét cellát, azaz 3.5 m-t fed le, így a rácson 3.6%-kal több víz jelenik meg a névleges terhelésnél. A geometriai változatban a cellákra jutó fedettséget a névleges terheléshez igazítottam. A zöld bozót 4.65 m-es sávja kilenc cellát (4.5 m) fed le.

**A hatásfok átskálázása.** A szükségleti tartományt $C^*$-hoz viszonyítva adtam meg, a lángok elérése pedig a szükséges sávszélességgel arányos ($W = 2.5\,L_f$). Átskálázáskor az $x/W$ és $y/W$ koordinátákon vett $C/C^*$ eloszlás mellett az $M/(W^2 C^*)$ arányt is meg kell őrizni. Ha a vízszintes hosszméreteket és $W$-t $b$-szeresére, a fedettséget és a szükségletet $a$-szorosára változtatom, ugyanazt a hatásfokot $ab^2$-szeres kibocsátott tömeggel kapom; a rögzített cellaméret és simítási hossz miatt a számított értékek kissé különbözhetnek.
Ha $b = 1$, az $aC$ fedettségű, $a$-szoros tömegű szóráskép hatásfoka $C^*$-nál megegyezik az eredeti szóráskép $C^*/a$-nál kapott hatásfokával, vagyis a hatásfokgörbe csak eltolódik a szükséglet tengelye mentén.

**A simítás és a cellarács (2.3. ábra).** Egy oszlop nedves szakasza mindig egész számú cellából áll. Ha 0.5 m-es cellákon egy egyenletes sáv pontosan a szükséges szélességű, és terhelése az ideális 1.5-szerese, diszkrét hatásfoka 0.67 vagy nulla, attól függően, hová esnek a szélei a cellaközéppontokhoz képest. A folytonos hatásfok a B-nél és a C-nél szinte állandó, az A-nál kissé ingadozik, 0.75–1 m-nél a legnagyobb, de ingadozása a saját $s_B$ szórásán belül marad. Zöld bozótban az A és a B olyan közel van egymáshoz, hogy 0.25 m-nél a folytonos hatásfok szerinti sorrendjük is megfordul; ennyire közeli szórásképeket nem lehet rangsorolni.

# B. függelék: A fizikai modell részletei

## B.1 A folyadék szállítása a cellák között

A megoldó beállításait a C.3. táblázat foglalja össze [OpenCFD2025].

| Jel | Jelentés |
|---|---|
| $f$ | a $P$ tulajdonos (owner) és az $N$ szomszéd (neighbour) cella közös cellafelülete |
| $\mathbf S_f$ | a $P$-ből $N$ felé mutató felületvektor |
| $s_{Pf}$ | $+1$, ha a cellafelület tulajdonosa $P$, egyébként $-1$ |
| $\phi_f$ | $\mathbf u_f\cdot\mathbf S_f$ |
| $H$ | folyadékban 1, levegőben 0 |

**Térfogatmérleg.** A (3.2) egyenlet a $P$ cellára és egy időlépésre integrálva, a Gauss–Osztrogradszkij-tétellel, közelítés nélkül [Roenby2016, 1.5–1.7. egyenlet]:

$$
\alpha_P^{n+1}=\alpha_P^{n}-\frac{1}{V_P}\sum_f s_{Pf}\,\Delta V_f,\qquad \Delta V_f=\int_{t^n}^{t^{n+1}}\!\!\int_f H\,\mathbf u\cdot d\mathbf S\,dt. \tag{B.1}
$$

Minden $\Delta V_f$ az egyik cellából távozik, a másikba érkezik, így a folyadék pontatlan $\Delta V_f$ mellett is megmarad; az egyetlen közelítés a $\Delta V_f$ becslése.

**isoAdvector** (`interIsoFoam`) [Roenby2016, 2–3. szakasz]. A módszer minden határfelszínt tartalmazó cellában megkeres egy sokszöget (izolapot), amelynek folyadék felőli oldalán pontosan $\alpha_PV_P$ folyadék van, és az időlépés alatt a normálisa mentén mozgatja. A cellafelület folyadékba eső $A_f(\tau)$ területe így minden pillanatban ismert:

$$
\Delta V_f\approx\frac{\phi_f}{|\mathbf S_f|}\int_{t^n}^{t^{n+1}}A_f(\tau)\,d\tau. \tag{B.2}
$$

A (B.2) közelítés egyenletes, időben állandó sebességtérben pontos [Roenby2016, 2.1–2.5. egyenlet].

**MULES** (`interFoam`, finom B futtatás). A MULES-séma a felszínt nem rekonstruálja, hanem mesterséges kompressziós fluxussal és korlátozással két-három cella szélesen tartja [Deshpande2012].

## B.2 A realizálható $k$–$\varepsilon$ turbulenciamodell

A Rouaix-eset cikkéhez [Rouaix2023] hasonlóan a realizálható $k$–$\varepsilon$ modellt [Shih1995] használtam, a keverékre egyetlen, a $\rho\mathbf u$ tömegfluxussal szállított $k$ és $\varepsilon$ mezővel.

**Transzportegyenletek.**

$$
\begin{aligned}
\frac{\partial(\rho k)}{\partial t}+\nabla\cdot(\rho\mathbf u k)&=\nabla\cdot\Big[\rho\Big(\nu+\frac{\nu_t}{\sigma_k}\Big)\nabla k\Big]+\rho\nu_t D^2-\rho\varepsilon,\\
\frac{\partial(\rho\varepsilon)}{\partial t}+\nabla\cdot(\rho\mathbf u\varepsilon)&=\nabla\cdot\Big[\rho\Big(\nu+\frac{\nu_t}{\sigma_\varepsilon}\Big)\nabla\varepsilon\Big]+\rho C_1 D\varepsilon-\rho C_2\frac{\varepsilon^2}{k+\sqrt{\nu\varepsilon}},
\end{aligned} \notag
$$

ahol $\nu_t=\mu_t/\rho$, $\mathbf D=\tfrac12[\nabla\mathbf u+(\nabla\mathbf u)^{\mathsf T}]$ és $D=\sqrt{2\mathbf D:\mathbf D}$. A jobb oldalon sorrendben a diffúziós, a keltési és a nyelőtag áll.

**Turbulens viszkozitás és együtthatók.** $\mu_t=\rho C_\mu k^2/\varepsilon$, ahol $C_\mu$ az alakváltozási sebességgel csökken (innen a „realizálható” jelző):

$$
C_\mu=\frac{1}{A_0+A_s U^{*}k/\varepsilon},\qquad C_1=\max\Big[0.43,\frac{\eta}{\eta+5}\Big],\quad \eta=\frac{Dk}{\varepsilon}, \tag{B.3}
$$

ahol $\mathbf W=\tfrac12[\nabla\mathbf u-(\nabla\mathbf u)^{\mathsf T}]$, $U^*=\sqrt{\mathbf D:\mathbf D+\mathbf W:\mathbf W}$, és $A_s$ a $\mathbf D$ invariánsaitól függ [Shih1995].

| Állandó | Érték | Forrás |
|---|---|---|
| $A_0$, $C_2$, $\sigma_k$, $\sigma_\varepsilon$ | 4, 1.9, 1, 1.2 | OpenCFD-alapértelmezés |

**Korlátok.** A modell a határfelszínen csak durva közelítés. A sűrűségugrás két oldalán ugyanaz a $k$ és $\varepsilon$ érvényes; a szabad felszín nem csillapítja a turbulenciát; az elkenődött határfelszínen $k$ keletkezhet; az $\alpha$ egyenletéből hiányzik a turbulens diffúziós tag.

## B.3 A cseppaprózódási modell

**Utódcsomagok.** Ha egy $n$ darab $d$ átmérőjű cseppet vivő cseppcsomag $K$ egyenlő tömegű utódcsomagra válik szét, a $d_c$ átmérőjű cseppeket vivő utódcsomag cseppszáma

$$
n_c=n\,(d/d_c)^3/K \notag
$$

így $\sum_c n_c d_c^3=nd^3$, a tömeg megmarad.

| Elem | Beállítás | Megjegyzés |
|---|---|---|
| Első aprózódás | méretosztályonként hat utódcsomag, $\mathbf u_p\pm\delta\,\mathbf e_j$ sebességgel | $\mathbf e_j$ véletlenszerűen elforgatott ortonormált bázis; a lökések összege nulla, az impulzus megmarad |
| Lökéssebesség | $\delta=s(\rho_g/\rho_l)^{1/2}=d/t^*$ | saját választás, mért töredéksebesség nem áll mögötte |
| További aprózódás | nincs újabb lökés | |
| Méretosztályok | öt Hsiang–Faeth-féle osztály, egy gyöknormális (root-normal) térfogateloszlás egyenlő tömegű tartományai; $\mathrm{MMD}=1.2\,\mathrm{SMD}$ | egyik sem nagyobb az eredeti cseppnél |

*Eltérések a forrásoktól.*

| Pont | Forrás | Saját felírás |
|---|---|---|
| Időzítés | azonnali szétbontás, amint $\mathrm{We}>12$, a keletkező cseppek azonnali újravizsgálatával [Shi2026] | aprózódás $T_b$ idő alatt (3.10) |
| Sebességlökés | véletlen, az eredeti sebességre merőleges; az impulzus nem marad meg [Shi2026, 11–13. egyenlet, 1. algoritmus] | hat utódcsomag, nulla összegű lökéssel |
| $t^*$ | a keletkező csepp átmérőjével [Shi2026, 14. egyenlet] | az eredeti csepp átmérőjével |
| Viszkózus tag | $\nabla\cdot(\mu\nabla\mathbf u)$, hiányzik a $(\nabla\mathbf u)\cdot\nabla\mu$ tag [Shi2026, 3. egyenlet] | az OpenFOAM `divDevRhoReff` operátora tartalmazza |
| Együttható | 0.776 [Shi2026, 15. egyenlet] | 0.766 (3.9); ezzel $T^*$ a módok határán 0.2%-on belül folytonos |
| Turbulencia | SST $k$–$\omega$ [Shi2026, 6–7. egyenlet] | realizálható $k$–$\varepsilon$ (B.2. függelék) |
| Méreteloszlás | zsákos és többmódusú aprózódásra | nyírásos aprózódásra is |
| $d_s$ stabil méret | a lassuló töredékfelhő pályája mentén [Pilch1987] | pillanatnyi felső korlát |
| Relatív sebesség (3.11) | a kezdeti [Hsiang1992] | az aprózódás végén érvényes; rekurzívan a keletkező cseppekre is |

## B.4 A szél, a turbulencia és a pályaszámítás

**Szélprofil.** Logaritmikus a talaj feletti $z_g$ magasság szerint ($z_0=0.03$ m feltételezett érdességi magasság):

$$
u_w(z_g)=u_{\mathrm{ref}}\,\frac{\ln(z_g/z_0)}{\ln(z_{\mathrm{ref}}/z_0)}. \tag{B.4}
$$

A (3.6) egyenletben $\mathbf w$ a kibocsátási magasságban fújó szél.

**Turbulens ingadozás.** Minden cseppcsomag saját, Ornstein–Uhlenbeck-folyamatot követő $\mathbf u'$ ingadozó légsebességet kap, amelyet az utódcsomagok örökölnek ($\zeta_j$ standard normális véletlen szám):

$$
u'_j(t+\Delta t)=f\,u'_j(t)+\sigma_j\sqrt{1-f^2}\,\zeta_j,\qquad f=e^{-\Delta t/T_L}. \tag{B.5}
$$

| Mennyiség | Érték | Forrás / megjegyzés |
|---|---|---|
| Szórások (szélirány, függőleges, keresztirány) | $(\sigma_1,\sigma_2,\sigma_3)=(2.4,1.25,1.9)\,u_*$, $u_*=0.4\,u_{\mathrm{ref}}/\ln(z_{\mathrm{ref}}/z_0)$ | [Panofsky1984] |
| Lagrange-féle időskála | $T_L=0.5\,z_g/\sigma_2$; 42 m-en körülbelül 100 s (az osztóval mintegy 70 s) | [Hanna1982], a Coriolis-paramétert ($f$) tartalmazó $1+15fz_g/u_*$ osztó nélkül; a kód lokálisan számítja, a talaj közelében rövidül (1 m-en körülbelül 2.4 s) |
| Függőleges diffúziós együttható | $K_2=\sigma_2^2T_L=0.625\,u_*z_g$ | a semleges felszíni réteg $0.4\,u_*z_g$ értékének mintegy 1.6-szorosa |
| Pályakeresztezés [Csanady1963] | nincs korrigálva | a csepp esés közben körülbelül két $0.5\,z_g$ méretű örvényen halad át; becslésem szerint az elmozdulás szórása 25–40%-kal kisebb volna, így (B.5) túlbecsli a szétszóródást |

**Légköri rétegződés.** A jelentés [Solarz2001, 5. ábra] csak a 2.24 m/s-os szelet közli (6.1 m-en [Suter2000, 9. o.]), ezért a 4.2. szakasz cseppcsomag-számításait feltételezett légkörökben is lefuttattam, szélirányonként egy véletlen realizációval, a 60%-ig töltött tartály semleges és labilis esetét pedig tíz turbulencia-realizációval is. A szél- és a turbulenciaprofilt a Monin–Obukhov-elmélet Businger–Dyer-féle alakjából vettem, $T_L$ mindenhol semleges alakú.

| Eset | Beállítás | Eltérés a semlegestől (hossz, 1–2 GPC-s terület) |
|---|---|---|
| Semleges (indifferens) | (B.4)–(B.5); a szél 42 m-en 3.05 m/s | – |
| Mérsékelten labilis | $L_{MO} = -30$ m, 1000 m vastag keveredési réteg; a $\partial K_z/\partial z$-hez hasonló sodródási tag elhanyagolva (becslésem szerint 1–2 cm/s) | 6%-on belül; 6 GPC-s terület cellánként egy cseppcsomaggal (60%-ig töltött tartály, tíz realizáció) 48.8 ± 1.4 m² helyett 20.2 ± 1.8 m² |
| Mérsékelten stabil | $L_{MO} = +60$ m; a szél 42 m-en körülbelül 4.1 m/s | 4%-on belül; a szórásképek 1–2 m-rel szélesebbek |
| Kanyargó szél | semleges; szélirány 20° szórással, 60 s korrelációs idővel | 6.5%-on belül |

Az öt szélirány átlagát véve a tíz realizáció között a hossz szórása 0.22–0.45 m, a szélességé 0.07–0.13 m; az 1 és 2 GPC-s területek relatív szórása cellánként egy cseppcsomaggal legfeljebb 1%, darabonként egy gömbbel legfeljebb 2.4%. Labilis légkörben a 6 GPC-s terület mind a tíz realizációban csökkent, cellánként egy cseppcsomaggal 24–31 m²-rel, darabonként egy gömbbel 3–14 m²-rel (79.6 ± 2.6 m² helyett 70.1 ± 3.5 m²). Az átlagos szélesség cellánként egy cseppcsomaggal minden realizációban 0.6–0.8 m-rel nőtt, darabonként egy gömbbel átlagosan 0.02 m-rel.

**Időintegrálás.**

| Lépés | Módszer |
|---|---|
| Sebesség, helyzet | explicit; előbb a sebesség, majd az új sebességgel a helyzet (párolgással félimplicit közegellenállás) |
| Időlépés | párolgás nélkül a legrövidebb közegellenállási relaxációs idő tizede; párolgással a relatív sebesség változási időskálájának tizede (a 0.002 s-nál rövidebb relaxációs idejű cseppek nélkül); mindkét esetben $10^{-5}$ és $5\times10^{-3}$ s között |
| Aprózódás előrehaladása (3.10) | a lépés elején kiértékelt $\Delta t/T_b$-vel nő, ha a lépés elején és végén is $\mathrm{We}>12$ |
| Aprózódás és talajra érkezés | ellenőrzés a lépés végén |

## B.5 A cseppek párolgása

Az izolált cseppre felírt modell [Denner2026, 12–33. egyenlet] az $r=d/2$ cseppsugár és a $T_d$ csepphőmérséklet változását adja meg:

$$
\frac{dr}{dt}=\frac{\mathrm{Sh}\,D_v M_w}{2r\rho_l\mathcal R}\left(\frac{\mathrm{RH}\,p_s(T_a)}{T_a}-\frac{p_s(T_d)}{T_d}\right),\qquad
\frac{dT_d}{dt}=\frac{3}{r\,c_l}\left[\frac{\mathrm{Nu}\,k_a\,(T_a-T_d)}{2r\rho_l}+L_v\frac{dr}{dt}\right]. \tag{B.6}
$$

Itt $T_a$ és $\mathrm{RH}$ a levegő hőmérséklete és relatív páratartalma; $p_s$ a telítési gőznyomás; $D_v$ a vízgőz diffúziós együtthatója; $k_a$ a levegő hővezetési tényezője; $M_w$, $\mathcal R$, $L_v$ és $c_l$ rendre a moláris tömeg, a gázállandó, a párolgáshő és a fajhő.

**Sherwood-szám.** A kód [Denner2026, 21. egyenlet] szerint számítja a csepp $\mathrm{Re}$ Reynolds-számával (3.1. táblázat) és $\mathrm{Sc}=\nu_g/D_v$-vel:

$$
\mathrm{Sh} = 2 + 0.316\,\mathrm{Re}\,\mathrm{Sc}^{2/3}\ \ \text{ha}\ \mathrm{Re}^{1/2}\mathrm{Sc}^{1/3} < 1.4,\qquad \mathrm{Sh} = 1.56 + 0.616\,\mathrm{Re}^{1/2}\mathrm{Sc}^{1/3}\ \ \text{egyébként}. \notag
$$

A két ág 0.216-os együtthatóval csatlakozna folytonosan, a közölt 0.316 valószínűleg sajtóhiba. A közölt értéket használtam; az alsó ágban az így kapott Sherwood-szám legfeljebb mintegy 8%-kal nagyobb, mint 0.216-os együtthatóval volna.

**Nusselt-szám.** Ugyanígy számítható, $\mathrm{Sc}$ helyén a levegő Prandtl-számával. A modell a nedves hőmérő hőmérsékletén alapuló korrekcióval egészíti ki [Denner2026, 22–25. egyenlet], és a kód ezt is alkalmazza.

# C. függelék: A szimulációk adatai és beállításai

*C.1. táblázat.* A dolgozatban felhasznált szimulációk. A VOF-futtatások számítási költsége a C.3. táblázat utolsó sorában szerepel.

| Vizsgálat | Hol szerepel | Megoldó vagy módszer | Fő beállítások |
|-----|---|-----|------------|
| Rouaix, finom A és B | 4.1. szakasz | OpenFOAM v2512, `interIsoFoam` (A), `interFoam` (B) | 10.1 millió cella, a háló sűrített tartománya a sugár pályája mentén 14 mm-es; C.3. táblázat |
| Rouaix, durva háló és változatai | 4.1. szakasz, 4.1. táblázat | OpenFOAM v2512, `interIsoFoam` | 2.73 millió cella, 0–5 s és 0–1 s; a változatok egy-egy beállításban térnek el |
| Black Hawk–Bambi (két peremfeltétel-változat és forgószárnyas) | 4.2. szakasz | `interIsoFoam`, majd a 3.6. szakasz cseppcsomag-modellje | feltételezett Torricelli-féle leürülés; teli tartály $C_d = 0.6/0.8/1.0$ (elsődleges 0.8), 60%-ig töltött tartály; 4.93 millió cella; rotoros futtatások impulzustárcsával, tartálytesttel, nagyobb tartományban |
| Sebességillesztett vízkibocsátások | 5.3–5.6. és 5.7. szakasz | mint a Bambi-futtatás | sugarak: négy kifolyónyílás, 2.498 m³ 3.5 s alatt, 15.0 millió cella; függöny: egy kifolyónyílás, 10.66 m³ 3.5 s alatt, körülbelül 5.2 millió cella |
| Sebességillesztett tervezési sávok | 5.6. szakasz, 5.2. táblázat | a 3.6. szakasz cseppcsomag-modellje | cseppméretenként 30 000 cseppcsomag (60 mm-nél 1500 kiinduló), $d_{50}$ = 1, 2, 4, 60 mm, 42 m, 42 kt, egyenletes és turbulens szél, öt szélirány |
| Az I4F redukált hálózati modellje | csak ebben a táblázatban | közönséges differenciálegyenletek, implicit középponti módszer | koncentrált paraméterek (a csőágak folyadékának tehetetlensége, a gázpárna nyomása); a szelepveszteség hasonló szerelvények alapján felvett változatai |

## A VOF-számítások beállításai

*C.2. táblázat.* A VOF-futtatások peremfeltételei. A légi járműhöz rögzített koordináta-rendszer origója a fúvóka, illetve a kifolyónyílás középpontja; $x$ a beáramló levegő irányába, $y$ felfelé mutat; a nehézségi gyorsulás $-9.81\,\mathbf e_y$ m/s².

| Perem | Sebesség | Nyomás ($p_{rgh}$) | $\alpha$ | $k$, $\varepsilon$ |
|-----|------|-----|--|-----|
| A levegő belépő pereme | rögzített: 70 m/s (Rouaix), 21.607 m/s (42 kt, a többi futtatás) | `fixedFluxPressure` | 0 | rögzített, az 1%-os intenzitásból |
| Nyitott peremek: Rouaix és a Bambi első változata | `pressureInletOutletVelocity` | `prghTotalPressure`, $p_0 = 0$, a dinamikus tag csak beáramlásnál | `inletOutlet` | `inletOutlet` |
| Nyitott peremek: a Bambi végleges változata, a forgószárnyas és a sebességillesztett futtatások | `freestreamVelocity` | `freestreamPressure` | `inletOutlet` | `inletOutlet` |
| A légi jármű fala | tapadási feltétel | `fixedFluxPressure` | 90°-os peremszög | falfüggvények, prizmatikus réteg nélkül (a Rouaix-fúvóka közelében $y^+ \approx 4000$) |
| Fúvóka vagy kifolyónyílás | rögzített: Rouaix 10 m/s lefelé; tartály: a leürülés időbeli táblázata; sebességillesztett: $(U, -v_n, 0)$ | `fixedFluxPressure` | 1 | rögzített értékek a Rouaix-fúvókán az 1%-os intenzitásból, a többi esetben a Bambi-ürítés átlagsebességéből |

| Mennyiség | Érték | Forrás / megjegyzés |
|---|---|---|
| Turbulenciaintenzitás, a Rouaix-fúvókán belépő víz | 1% | [Rouaix2023, 3. táblázat]; a többi futtatás kifolyónyílásán a Bambi-ürítés 2.74 m/s-os átlagsebességéből számított, rögzített $k$ és $\varepsilon$ |
| Turbulenciaintenzitás, belépő levegő | 1% | feltevés |
| Turbulens hosszlépték | $0.07\,d$ | feltevés; 70 m/s-nál $k = 0.735$ m²/s², $\nu_t/\nu \approx 830$ |

**A nyitott peremek két változata (tartályos eset).**

| Mennyiség | Érték | Megjegyzés |
|---|---|---|
| Hidrosztatikus nyomáskülönbség, 32 m magas tartomány | 370 Pa | nagyobb a dinamikus nyomásnál |
| Dinamikus nyomás | 275 Pa | |
| Első változat (össznyomás-perem) | a levegőoszlop lesüllyedt és alul kiáramlott | ezzel készültek a $C_d = 0.6$ és $1.0$ futtatások |
| Végleges változat (szabadáramlási perem) | nincs lesüllyedés | 20 szórásképen átlagosan hossz +0.6%, szélesség −1.0%; átalakított víz 63.2% helyett 64.0% |

*C.3. táblázat.* A VOF-futtatások beállításai. Célérték: beállított; elért: a futtatási naplóból (finom B: a megőrzött beállítófájlokból).

| | Rouaix 1. eset (4.1. szakasz) | Black Hawk–Bambi (4.2. szakasz) | Sebességillesztett kibocsátások (5. fejezet) |
|---|-----------|------|------|
| Megoldó | OpenFOAM v2512, `interIsoFoam` (isoAdvector); finom B: `interFoam` (MULES) | `interIsoFoam` | mint a Bambi |
| Közegek | víz 997.6 kg/m³, levegő 1.18 kg/m³, $\sigma = 0.072$ N/m [Rouaix2023, 2. táblázat] | mint a Rouaix | mint a Rouaix |
| Turbulencia | realizálható $k$–$\varepsilon$ (B.2. függelék); a durva hálón lamináris változat is | realizálható $k$–$\varepsilon$ | mint a Bambi |
| Számítási tartomány | $x$ irányban 20 m, mélysége 10 m, szélessége felül 7 m, alul 16 m; a belépő perem 2.5 m-rel a fúvóka előtt; a fúvóka körül feltételezett, 0.08 m-es Gauss-alakú falbemélyedés | $x$ = −3.46 … 44.93 m, $y$ = −32.26 … 0 m, $z$ = ±14.4 m; rotoros futtatások: nagyobb tartomány, feltételezett hengeres tartálytest | mint a Bambi; $S$ 6 m-rel a kifolyónyílások alatt ($x < 28$ m, $y > -6$ m) |
| Kibocsátás | $d_j = 0.40$ m, egyenletes, 10 m/s lefelé | $d = 0.508$ m-es szelep, Torricelli-féle leürülés; elsődleges: teli tartály, $C_d = 0.8$ | $U$ sebességű, $h = a v_n/U$ vastag réteg ($a$ a kifolyónyílás hossza a nyomvonal mentén), feltételezett 0.1 s-os fel- és lefutással. Sugarak: négy 0.108 m széles kifolyónyílás 0.54 m-es osztással, 79 mm-es réteg. Függöny: egy 1.008 m széles kifolyónyílás, 0.144 m-es réteg |
| Levegő | 70 m/s | 21.607 m/s (42 kt); a VOF-számításban nincs szél | mint a Bambi |
| Dimenziótlan számok | $q = 17.3$, $\mathrm{We}_g = 3.2 \times 10^4$, $\mathrm{Fr} = 5.0$ | $q = 13.6$, $\mathrm{We}_g = 3.9 \times 10^3$, $\mathrm{Fr} = 1.2$ (átlagos kilépési sebesség 2.74 m/s) | a nyomvonal mentén nincs relatív sebesség |
| Háló | durva: 2.73 millió cella, tenzorszorzat-háló, a sugár körül $\Delta x$, $\Delta y$ 34.5–69 mm, $\Delta z$ 44–67 mm. Finom: 10 099 265 cella, a közölt pálya mentén $y \approx -2.5$ m-ig 14 mm, azon túl 28–112 mm; nemortogonalitás legfeljebb 58° | 4.93 millió hexaéder, lépcsőzetes (castellated) háló: a kifolyónyílásnál 18 mm, $y = -5$ m-ig 36 mm, $\{x < 5, y > -10\}$ m-ben 72 mm, azon túl 144–576 mm; nemortogonalitás legfeljebb 25° | 18, 36, $S$ felé 72, azon túl 144–576 mm. Sugarak: 15.0 millió cella, $S$-ig 36 mm-es sűrített tartománnyal (a sugár 18, 36 és 72 mm-en 4.4, 2.2, illetve 1.1 cella vastag). Függöny: körülbelül 5.2 millió cella (8, 4, illetve 2 cella). Nemortogonalitás körülbelül 25° |
| Sémák | Euler; gradiens `Gauss linear`; $\mathbf u$: `limitedLinearV 1`; $k$, $\varepsilon$: upwind; korrigált Laplace-tag. Finom B: $\alpha$: `vanLeer`, kompresszió `linear`, $c_\alpha = 1$ | mint a Rouaix | mint a Bambi |
| Határfelszín | `isoAlpha`, `isoFaceTol` és `surfCellTol` $10^{-10}$, `nAlphaBounds` 5, `clip` be (nem térfogatmegőrző), `snapTol` 0, egy alciklus. Finom B: `MULESCorr`, 2 korrektor, 3 alciklus, 5 korlátozó-iteráció | mint a Rouaix | mint a Bambi |
| Lineáris megoldók | $p_{rgh}$: PCG + GAMG, $10^{-8}$, relatív 0.01 (utolsó korrekciónál 0); $\mathbf u$: PBiCGStab/DILU, $10^{-8}$, relatív 0.1. Finom A: korrekciónként 1–2 iteráció, az utolsónál hat | $p_{rgh}$: PCG + DIC, azonos toleranciák | mint a Bambi |
| PIMPLE (külső iteráció, nyomáskorrekció, nemortogonális korrekció) | 1, 3, 1. Finom B: 3, 2, 1, relaxáció 0.7 ($p_{rgh}$), 0.8 ($\mathbf u$, $k$, $\varepsilon$) | 1, 3, 1 | 1, 3, 1 |
| Időlépés, célérték (legnagyobb $\mathrm{Co}$, $\mathrm{Co}_\alpha$, $\Delta t$) | 1, 0.5, 0.2 ms. Finom B: 5, 2, 1 ms | 1, 0.5, 1 ms | 1, 0.5 |
| Időlépés, elért érték | finom A, 1092 lépés 1.42 és 1.50 s között: legnagyobb $\mathrm{Co}$ átlaga 0.85, csúcsa 1.10; $\mathrm{Co}_\alpha$ 0.48, illetve 0.80; átlagos időlépés 73 μs. A célértékek felezése a durva hálón legfeljebb 2.2%-kal változtatta a behatolási mélységet | | |
| Kimenet | 5 ms-onként a sugár mélysége és szélessége vékony szeletekben; folyadékmérleg | minden időlépésben a folyadéktérfogat és a pereken, valamint $S$-en átáramló térfogat; $S$-en túli cellasávok 0.1 s-onként; 0–8 s | ugyanazok a kimenetek, 0–5 s |
| Folyamatok, költség | durva: 120 folyamat, 5 s-ra 619 processzormag-óra. Finom A: 0.4248 s után 240 folyamat, szimulált másodpercenként körülbelül 1850 processzormag-óra | 120 folyamat, 8 s-ra 1219 processzormag-óra | 240 folyamat, a két futtatás együtt körülbelül 4600 processzormag-óra |

**A Bambi-tartály kibocsátása.**

| Mennyiség | Érték | Forrás / megjegyzés |
|---|---|---|
| Repülési sebesség | 42 kt (21.607 m/s) | mérési jelentés [Solarz2001, 1. o. és 5. ábra] |
| A tartály alja a talaj fölött | 42.06 m | ugyanott |
| Szél | 2.24 m/s, 6.1 m-en | ugyanott; mérési magasság [Suter2000, 9. o.] |
| Víz a tartályban | legfeljebb 2.50 m³, 53 hüvelyk (1.35 m) töltési magasságig | ugyanott |
| Szelep | 20 hüvelykes (0.508 m-es), „egyetlen térfogatárammal” | ugyanott |
| Kibocsátott vízmennyiség | – | nincs közölve |
| Tartály és ürülés | prizmatikus tartály, Torricelli-féle ürülés | feltevés |

Az $A_v$ szelepkeresztmetszettel, $h_0$ kezdeti töltési magassággal és $V_0$ kezdeti térfogattal:

$$
Q(t) = Q_0(1 - t/T),\qquad Q_0 = C_d A_v \sqrt{2 g h_0},\qquad T = 2V_0/Q_0,
$$

az elsődleges futtatásban $Q_0 = 0.833$ m³/s és $T = 6.0$ s. A C.3. ábra tartályváltozataiban a 2.498 m³, az 1.346 m-es töltési magasság és a 0.508 m-es szelep változatlan, csak a tartály bővül felfelé.

*(Ábra: `bambi_setup.png`, lásd a dolgozat PDF-változatát.)*
*C.1. ábra.* A Bambi-tartály kibocsátásának szimulációja (teli tartály, $C_d = 0.8$, 42 kt): (a) oldal- és (b) elölnézet a tartománnyal és a peremfeltételekkel; (c) előírt és kilépő térfogatáram; (d) cellaméret a középsíkon a $\chi = 0.5$ szintvonallal és az $S$ átalakítási felülettel; (e) a kifolyónyílás környezetének cellái.

## A Rouaix-eset mélységgörbéiből számított gyorsulás

Az alábbi táblázat mélységgörbénként megadja a görbét létrehozó állandó áramlásirányú gyorsulást.

| Görbe | Gyorsulás [m/s²] | Megjegyzés |
|---|---|---|
| Saját görbék | 25–32 | a $\chi = 0.1$ felhő itt 2.5–4 sugárátmérő széles (4.3b ábra); ekkora szélességnél a közegellenállás szokásos ellenállástényezővel is ennyit ad |
| A [Rouaix2023] cikk 13. ábrája, $x$ a fúvóka áramlás felőli szélétől | 48–59 | |
| Ugyanaz, $x$ a fúvóka középpontjától | akár 176 | a teljes sugárra ható közegellenállásból adódó gyorsulás 16-szorosa |

## Az átalakítás részletei

**Darabonként egy gömb.** Az egymással érintkező áthaladó cellák egy folyadékrészt alkotnak, ebből egyetlen cseppcsomag lesz.

| Mennyiség | Meghatározás |
|---|---|
| Térfogat | a cellák $\alpha_cV_c$ folyadéktérfogatainak összege |
| Helyzet, sebesség | a cellák helyzetének és sebességének e térfogatokkal súlyozott átlaga |
| Csepp | egyetlen, $d_{eq}$ átmérője az azonos térfogatú gömbé |

A cseppcsomag pontosan a cellák folyadékát és impulzusát viszi tovább. A (3.5) feltétel a cella teljes folyadéktérfogatát egyben számolja el, így pillanatfelvételből csak becsüli az $S$-en áthaladó folyadékot.

**A $k_i$ szorzó.** Megőrzi a folyadék térfogatát minden olyan mentési időközben, amelyben legalább egy cseppcsomag keletkezik. Az átmérőket $k_i^{1/3}$-szorosukra növeli (a tartályos kibocsátásban mintegy negyedével), a helyzetet és a sebességet nem változtatja. Nem őrzi meg az $S$-en áthaladó impulzusáramot, az időközön belüli időzítést és a kihagyott folyadékrészek méretét.

**A folyadékrészek mérete az $S$-en (4.6e ábra).** 813 folyadékrész, skálázás előtt, $\alpha \ge 0.05$; a 0.04–0.1, 0.1–1 és 1–10 m-es osztályba rendre 340, 462 és 0 darab esik, további 11 darab 0.04 m-nél kisebb (ugyanezekkel az osztályokkal jellemezte a folyadék 1 s-os állapotát a CL-415 és a Dash-8 kibocsátásának számítása [Calbrix2023, 11. ábra]). A darabszám szerinti medián átmérő 0.107 m (1.5 cella), a térfogat szerinti 0.20 m (2.8 cella), és a térfogat 93%-a a 0.1–1 m-es osztályba esik. Az oldalirányú szétterülési sebesség mindkét osztályban körülbelül 0.9 m/s, ami nagyságrendben egyezik a [Calbrix2023] eredményével.

**Vízmérleg az $S$-en (teli tartály).** A 4.6d–f ábra időpontjaiban (4, 9, 20 s) a 2488 kg vízből rendre 1, 1039 és 2254 kg ért talajt. Minden vizet tartalmazó cellára alkalmazva a (3.5) szabály a víz 64.0%-át alakítja át (hat tartályos futtatás: 60–69%, annál kevesebb, minél kisebb a térfogatáram). Egy 0.3 s hosszú, 0.01 s-onként mentett újraindításban az $S$-nél a víz 84–92%-a félig sem kitöltött ($\alpha<0.5$) cellákban van, és a megoldó 1.34-szer annyi vizet szállít át, mint amennyi a cellafelületi térfogatáram és az áramlás felőli cella fázisarányának szorzatából adódna; a szabály a 0.01 s-onként mentett mezőkből a víz 75%-át, a 0.1 s-onként mentettekből 61%-át adja vissza. A víz 0.02–0.2%-a olyan mentési időközökben halad át, amelyekben nem keletkezik cseppcsomag; ezt nem alakítottam át (a 4.7. ábra futtatásában 0.2%, ezért ott a cseppcsomagok 2494 helyett 2488 kg vizet visznek).

## A Bambi-szórásképek kiértékelése

*C.4. táblázat.* A 4. fejezetben használt mennyiségek meghatározása.

| Mennyiség | Meghatározás |
|------|------------------------|
| Rouaix-féle behatolási mélység, szélesség | a legmélyebb $\alpha \ge 0.1$ cella, illetve a keresztirányú kiterjedés a mérési helyek körüli 6 cm-es szeletekben, 5 ms-onként; átlag 0.68–1.50 s (finom A, B) vagy 0.6–1.0 s (durva). A durva hálóról közölt egyes mélység- és szélességértékek vékony szeletek helyett síkmetszetekből származnak. A cikk [Rouaix2023] 12a ábrája a fúvókakörhöz viszonyítva leolvasva (±3%, 1σ) |
| Folyadékmérleg; átalakított hányad | bejuttatott = tartományban lévő + kiáramlott (naplózott értékek); a (3.5) szabállyal cseppcsomaggá alakított térfogat osztva az $S$-en átáramló naplózott térfogattal |
| Folyadékrészek $S$-nél | 4.6e ábra: végleges változat, $\alpha \ge 0.05$, skálázás előtt; 4.6a–d és 4.7. ábra: első változat (azonos háló és kibocsátás) |
| Gyűjtőedényes leolvasás | edények az $x_0 + k\Delta$, $y_0 + l\Delta$ pontokban, $\Delta$ = 3.048 m (10 láb), a 0.5 m-es cellájuk értékével; lineáris visszainterpolálás a cellákra, lefelé kerekítés a 0.5, 1, 2, 3, 4, 6, 8 és 10 GPC-s szintekre, szintközi interpoláció (2.2. szakasz); 16 rácshelyzet átlaga, $x_0, y_0 \in \{0, \Delta/4, \Delta/2, 3\Delta/4\}$ |
| Hossz, szélesség | a víztömeg 5–95%-át tartalmazó szakasz a legalább 0.5 GPC-t kapott részen; a 4.2. táblázatban szintközi, a C.5. táblázatban lépcsős leolvasással (a hosszak körülbelül 1%-kal térnek el) |
| Bemutatott rácshelyzet | 4.9. és C.4. ábra: a 16 helyzet átlagához legközelebbi (4.9. ábra: $x_0 = 0.76$ m, $y_0 = 2.29$ m) |
| Legjobb szóráskép | az egyenletes és a semleges légköri szélben számított 40 közül az, amelyre $\lvert\mathrm{számított}/\mathrm{mért} - 1\rvert$ összege a hosszra, a szélességre és az 1, 2, 6 GPC-s területre a legkisebb; a módot az összevetés előtt rögzítettem, a további 60 szóráskép csak ezután készült |
| Leolvasási tartomány | a jelentés [Solarz2001] 5. ábrájának összes leolvasását (szintközi, lépcsős, szűkebb és tágabb digitalizálás) tartalmazó tartomány, bővítve a 16 rácshelyzet átlagos fél terjedelmével, 0.5–2 GPC-n a 2. táblázat 84 kt-s és a 7. ábra 81 kt-s, eltérő magasságú Black Hawk-kibocsátásának eltérésével |
| Az ismételt kibocsátások tartománya | a leolvasási tartomány $e^{\pm c}$-szerese; $c$ a merevszárnyú kibocsátások és átlaguk logaritmikus eltéréseinek négyzetes középértéke ([McFayden2023] A2-6. ábra, 68–71. o., és A2-2. táblázat): a hosszra, a szélességre és a területet közelítő hossz–szélesség szorzatra 0.5–1 GPC-n 0.19, 0.14, 0.29; 2–3 GPC-n 0.68, 0.21, 0.80; 4 GPC-től 1.2, 1.2, 2.4; csak becslés (nem ismételt mérések, külső függesztésű tartály nélkül). A keresztirányú csúcs/átlag aránynál a három Black Hawk-kibocsátás [Solarz2001, 5–7. ábra] is benne van |
| Tájolás | C.6. ábra: a több vizet kapott szélső harmad jelzi a kezdetet; ez cellánként egy cseppcsomaggal mind a 12, darabonként egy gömbbel 7 szórásképen teljesül. A tömegközépponttól $\pm 0.6L$-en belül a fordított tájolás illeszkedik jobban a legjobb szórásképhez (négyzetes közepes eltérés 0.193, szemben 0.243-mal), a helyes a 20 egyenletes szélű szórásképhez (0.235, szemben 0.249-cel) |
| Leolvasási veszteség | 4.2. szakasz: a gyűjtőedényes leolvasással kapott és a 0.5 m szórású Gauss-függvénnyel simított számított vízeloszlás különbsége |

*C.5. táblázat.* Öt Bambi-változat összevetése a jelentés [Solarz2001] 5. ábrájával, 2.24 m/s-os egyenletes szélben; gyűjtőedényes, lépcsős leolvasás. Hossz, szélesség és talajszinti tömeg szerinti medián cseppátmérő ($d_{50}$): tíz szóráskép (öt szélirány, két cseppcsomag-felépítés) átlaga, zárójelben a terjedelem. $t_{50}$ (tömeg szerinti medián repülési idő): egyetlen oldalszeles, darabonként egy gömbbel számított pályasorozat. $U(t_{95}-t_{05})$: az $S$ felületen áthaladó víz középső 90%-ához tartozó idő és a repülési sebesség szorzata.

| Eset | $U(t_{95}-t_{05})$ [m] | Hossz [m] | Szélesség [m] | $t_{50}$ [s] | Talajszinti $d_{50}$ [mm] |
|------|--|--|--|--|--|
| Mért, 5. ábra | – | 65 | 10 | – | – |
| Teli tartály, előírt légáramlás | 101 | 99 (91–111) | 14.6 (6.8–24.2) | 1.94 | 1.16 |
| Teli tartály, saját légáramlás | 101 | 101 | 13.9 | 1.73 | 1.66 |
| Teli tartály, forgószárny, előírt légáramlás | 102 | 98 (86–113) | 14.8 (6.5–24.7) | 2.00 | 1.13 |
| Teli tartály, forgószárny, saját légáramlás | 102 | 98 (88–110) | 13.7 (7.0–21.5) | 1.71 | 1.58 |
| 60%-ig töltött tartály, előírt légáramlás | 77 | 75 (67–82) | 13.4 (6.7–21.4) | 2.92 | 0.87 |

$U(t_{95}-t_{05})$ a repülési sebesség és a folyadék középső 90%-a $S$-en való áthaladási idejének szorzata; $t_{50}$ és a talajszinti $d_{50}$ a tömeg szerinti medián repülési idő, illetve cseppátmérő a talajon; „saját légáramlás” a futtatás időátlagolt VOF-légsebessége.

**Rotoros futtatások.**

| Mennyiség | Érték | Megjegyzés |
|---|---|---|
| Modell | impulzustárcsa (actuator disk) | |
| Tolóerő | 97.9 kN | az UH-60A legnagyobb felszálló tömegéhez |
| Tárcsasugár | $R = 8.18$ m | |
| A tárcsa magassága a tartály fölött | 20.2 m | feltételezett 15 m-es (50 láb) kötél és 5 m |
| Indukált sebesség | 16–20%-kal nagyobb az impulzuselméleti értéknél | a nyomdok inkább túl erős; függésben négy, csak levegővel végzett próbafuttatással ellenőriztem, egyik sem teljesítette az előre rögzített ±10%-os feltételt |

A forgószárnyas eset eltérésében a tartálytest és a nagyobb tartomány hatása is benne van.

*(Ábra: `ch4_length_rule.png`, lásd a dolgozat PDF-változatát.)*
*C.2. ábra.* A Bambi-változatok gyűjtőedényes leolvasással kapott szóráskép-hossza az $U(t_{95}-t_{05})$ szorzat függvényében. A hibasávok futtatásonként húsz szóráskép terjedelmét jelzik (öt szélirány, két cseppcsomag-felépítés, egyenletes és semleges légköri szél). A kihúzott vonal az 1:1 arányt, a szaggatott a mért 65 m-t mutatja.

*(Ábra: `ch4_release_family.png`, lásd a dolgozat PDF-változatát.)*
*C.3. ábra.* A repülési sebesség és a víz középső 90%-a kifolyási idejének szorzata a tartály felső és alsó átmérőjének aránya függvényében, különböző alakú, töltöttségű és $C_d$ kifolyási tényezőjű tartályokra. Okkersárga vonal: a számított szóráskép hossza; szürke sáv: a mért, mintegy 65 m-es hosszt adó tartomány.

*C.6. táblázat.* A Black Hawkról Bambi 660-as tartállyal végzett három mért vízkibocsátás [Solarz2001, 5–7. ábra], a 4.2. táblázattal azonos leolvasással (ahol a lépcsős és az interpolált eltér, mindkettő szerepel). A hossz és a szélesség a talajra jutott víz 5–95%-át fogja közre; a vízmennyiséget és a szélirányt nem közölték.

| Mérés | Sebesség, $U$ [m/s] | Hossz, $L$ [m] | $L/U$ [s] | $3\ \mathrm{s} \times U$ [m] | Szélesség [m] | A szórásképen méterenként lévő víz [kg/m] |
|---|---|---|---|---|---|---|
| 5. ábra, 42 kt | 21.6 | 65–66 | 3.0–3.1 | 65 | 10 | 23–27 |
| 6. ábra, 50 kt | 25.7 | 88–89 | 3.4–3.5 | 77 | 14 | 17–21 |
| 7. ábra, 81 kt | 41.7 | 100 | 2.4 | 125 | 15 | 9–12 |


*(Ábra: `ch4_map_comparison.png`, lásd a dolgozat PDF-változatát.)*
*C.4. ábra.* A Black Hawk-kibocsátás mért szórásképe [Solarz2001, 5. ábra] (felül) és a teli, illetve 60%-ig töltött tartállyal számított szórásképek (cellánként egy cseppcsomag, gyűjtőedényes leolvasás) öt irányú egyenletes szélben és légköri szembeszélben. A számított szórásképek egy rácshelyzetet mutatnak, méreteik 16 rácshelyzet átlagai; pontozott téglalap: 65 × 10 m.

*(Ábra: `rotor_wake_centre_plane.png`, lásd a dolgozat PDF-változatát.)*
*C.5. ábra.* A levegő–víz keverék időátlagolt függőleges sebessége ($\bar u_y$) a forgószárnyas futtatásban (teli tartály, 42 kt). (a) A középsíkhoz legközelebbi mentett sík ($z = 0.23$ m); vastag vonal: a forgószárny tárcsája, pont: a szelep, szaggatott vonal: a $B$ tartomány határa. (b) Átlag az $5 \le x \le 45$ m, $|z| \le 3$ m tartományban, ahol a cseppcsomagok hullanak, a forgószárnyas futtatásban. A szürke függőleges vonal a nulla sebességet jelöli.

*(Ábra: `ch4_profiles.png`, lásd a dolgozat PDF-változatát.)*
*C.6. ábra.* A víz hossz- (balra) és keresztirányú (jobbra) eloszlása a Black Hawk-szórásképeken, a tömegközéppont körül a hosszal, illetve a szélességgel normálva; a kibocsátás balról indul. Szürke: a három vízzel végzett Black Hawk-mérés tartománya [Solarz2001, 5–7. ábra]; vastag okkersárga: a 4.9b ábra számított szórásképe. Felül: cellánként egy cseppcsomag; alul: darabonként egy gömb.

**Tartályváltozatok és a sebességsorozat.** A mért hosszt adja például a felfelé 1.3-szeresére bővülő, 80%-ig töltött tartály $C_d = 1.0$ mellett, vagy az 1.6-szeresére bővülő, 60%-ig töltött tartály $C_d = 0.8$ mellett (C.3. ábra). A nyomvonal egy méterére jutó víz 42-ről 81 kt-ra mintegy 25-ről 10 kg/m-re csökken (C.6. táblázat); azonos térfogatáramnál csak a felére csökkenne (5.2. szakasz). A rövidebb ürítési idő ezt önmagában nem magyarázza, ehhez az átlagos térfogatáramot és a kimutatott víz hányadát is ismerni kellene. A gyorsabb kibocsátást (80%-ig töltött tartály, $C_d = 1.0$) egyetlen alkalommal választottam ki, a 42 kt-s hosszhoz illesztve. A vele mindhárom sebességre, a teli tartállyal pedig 50 és 81 kt-ra elindított öt futtatás közül egyik sem jutott el a tervezett végidőig.

## Az 5. fejezet kiegészítő ábrái

*(Ábra: `ch5_credit.png`, lásd a dolgozat PDF-változatát.)*
*C.7. ábra.* A teli külső függesztésű tartály és a négy sebességillesztett sugár vízeloszlása sűrű bozótra értékelve (6 GPC, 3.38 m sávszélesség), szembeszélben, cellánként egy cseppcsomaggal. Jelölések, mint a 2.8. ábrán.

*(Ábra: `ch5_crosswind.png`, lásd a dolgozat PDF-változatát.)*
*C.8. ábra.* A teli tartály, a négy sugár és a vízfüggöny szórásképe 2.24, 5 és 10 m/s-os oldalszélben (42 kt, 42 m), a nyomvonalat tartva (szaggatott vonal), cellánként egy cseppcsomaggal. A felső panelek a szóráskép közepéből hosszirányban 76 m-es részletet mutatnak. Az alsó két grafikon a legalább 0.5 GPC fedettségű cellák vizét jellemzi. Balra a keresztirányú tömegeloszlás 5. és 95. percentilise közötti szélesség, jobbra a tömegközéppont nyomvonaltól mért eltolódása látható. A kihúzott vonal a cellánként egy cseppcsomaggal, a pontozott a darabonként egy gömbbel kapott eredményt jelöli.

*(Ábra: `ch5_video_vs_vof.png`, lásd a dolgozat PDF-változatát.)*
*C.9. ábra.* Az I4F eszköz habkibocsátásának képkockái (Mi-8, körülbelül 27 m/s, 11–13 m magasság), rajtuk az 5.4. szakasz számított vízsugara (kékeszöld; 21.6 m/s, 42 m, forgószárny nélkül), közelítő méterskálán. (d) A csóva mélysége a kifolyónyílás alatt. A videópontok hibasávjának félhossza a 10 képpontos leolvasási bizonytalanság méterben kifejezett értékének és a leolvasott mélység 10%-ának összege. (e) A csóva hátsó szélének lemaradása a nyílás mögött, $U t$-hez viszonyítva (az 1-es érték felel meg a sebességillesztett kibocsátásnak).
