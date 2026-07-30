# QLBM – quantum D2Q9 streaming prototype

Tento projekt ověřuje kvantovou implementaci **streamingového kroku D2Q9
Lattice Boltzmann Method (LBM)** pomocí Qiskitu. Podporuje ideální lokální
simulátor i skutečný IBM Quantum backend. Výsledky každého časového snímku
ukládá jako obrázek a jako numerická data pro následnou rekonstrukci hustoty,
rychlosti a tlaku.

Projekt obsahuje tři úrovně implementace:

1. původní basis-state experiment pro samostatné ověření `StreamG`;
2. klasický referenční D2Q9 BGK solver;
3. hybridní QLBM, ve kterém jsou collision a okrajové podmínky klasické,
   zatímco streaming úplného populačního pole probíhá kvantovým obvodem.

Hybridní varianta počítá density, velocity a pressure a podporuje ideální
Aer simulátor i experimentální IBM QPU backend. Nejde však o plně kvantový
CFD solver: BGK collision, inlet/outlet a makroskopická rekonstrukce stále
probíhají na klasickém procesoru.

## Matematický model D2Q9

LBM nepracuje přímo s Navierovými–Stokesovými rovnicemi. V každé buňce
\(\mathbf{x}=(x,y)\) uchovává devět distribučních funkcí

\[
f_i(\mathbf{x},t), \qquad i=0,\ldots,8,
\]

které reprezentují populaci částic pohybujících se diskrétními rychlostmi
\(\mathbf{c}_i\).

Použité směry jsou:

```text
6  2  5
 \ | /
3  0  1
 / | \
7  4  8
```

\[
\begin{aligned}
\mathbf{c}_0&=(0,0),&
\mathbf{c}_1&=(1,0),&
\mathbf{c}_2&=(0,1),\\
\mathbf{c}_3&=(-1,0),&
\mathbf{c}_4&=(0,-1),&
\mathbf{c}_5&=(1,1),\\
\mathbf{c}_6&=(-1,1),&
\mathbf{c}_7&=(-1,-1),&
\mathbf{c}_8&=(1,-1).
\end{aligned}
\]

V kódu je primární datový formát:

```python
f.shape == (9, ny, nx)       # f[direction, y, x]
```

### Makroskopické veličiny

Hustota je součet všech devíti populací:

\[
\rho(\mathbf{x},t)=\sum_{i=0}^{8} f_i(\mathbf{x},t).
\]

Hybnost:

\[
\rho\mathbf{u}
  =\sum_{i=0}^{8} f_i(\mathbf{x},t)\mathbf{c}_i.
\]

Rychlost pro \(\rho\neq0\):

\[
\mathbf{u}(\mathbf{x},t)
  =\frac{\sum_i f_i\mathbf{c}_i}{\rho}.
\]

V buňce s nulovou hustotou nastavuje implementace rychlost na nulu, aby
nevzniklo dělení nulou.

Pro izotermální D2Q9 platí rychlost zvuku v lattice jednotkách

\[
c_s^2=\frac{1}{3}
\]

a tlak se rekonstruuje jako

\[
p=c_s^2\rho=\frac{\rho}{3}.
\]

Výstupní rozměry:

```python
density.shape == (ny, nx)
velocity.shape == (ny, nx, 2)  # (..., ux, uy)
pressure.shape == (ny, nx)
```

Rekonstrukci implementuje
[`reconstruct_fields(f)`](src/macrocomputations.py).

### Kompletní klasický LBM krok

Plný BGK/SRT LBM časový krok by obsahoval collision a streaming:

\[
f_i^*(\mathbf{x},t)
=f_i(\mathbf{x},t)
-\frac{1}{\tau}
\left(f_i(\mathbf{x},t)-f_i^{eq}(\mathbf{x},t)\right),
\]

\[
f_i(\mathbf{x}+\mathbf{c}_i\Delta t,t+\Delta t)
=f_i^*(\mathbf{x},t).
\]

Rovnovážná distribuce D2Q9 je

\[
f_i^{eq}
=w_i\rho
\left[
1+\frac{\mathbf{c}_i\cdot\mathbf{u}}{c_s^2}
+\frac{(\mathbf{c}_i\cdot\mathbf{u})^2}{2c_s^4}
-\frac{\mathbf{u}\cdot\mathbf{u}}{2c_s^2}
\right],
\]

kde

\[
w_0=\frac{4}{9},\qquad
w_{1,2,3,4}=\frac{1}{9},\qquad
w_{5,6,7,8}=\frac{1}{36}.
\]

Kinematická viskozita v lattice jednotkách souvisí s relaxační dobou:

\[
\nu=c_s^2\left(\tau-\frac{1}{2}\right)\Delta t.
\]

Collision operator, \(f_i^{eq}\), viskozita a škálování lattice jednotek na
SI jednotky jsou zde uvedeny jako cílový fyzikální model, ale současný QPU
obvod je ještě neimplementuje.

## Implementovaný streaming

Projekt implementuje periodické zobrazení

\[
|x,y,d\rangle
\longmapsto
\left|
(x+c_{d,x})\bmod N_x,\,
(y+c_{d,y})\bmod N_y,\,
d
\right\rangle.
\]

Směr \(d\) se streamingem nemění. Periodické modulo znamená, že částice po
opuštění pravého/horního okraje vstoupí z opačné strany.

Obecná brána `StreamG` podporuje i superpozici směrů a používá direction
controls. Na dnešním QPU je však kvůli multi-controlled aritmetice velmi
hluboká.

Pro měřený basis-state experiment používá produkční cesta `BasisStreamG`.
Směr je před spuštěním známý, proto se vloží pouze potřebný increment nebo
decrement bez neaktivních direction větví. Jde o ekvivalentní streaming pro
známý basis stav, nikoliv o obecný superpoziční algoritmus.

## Kvantové kódování

Pro mřížku \(N_x\times N_y\), kde jsou rozměry mocninou dvou:

\[
n_x=\log_2N_x,\qquad n_y=\log_2N_y.
\]

Použité registry:

```text
|x> : REG_X_SIZE qubitů
|y> : REG_Y_SIZE qubitů
|d> : 4 qubity
```

Čtyři direction qubity kódují hodnoty 0–15, ale fyzikálně platné D2Q9 směry
jsou pouze 0–8. Výskyty 9–15 se proto reportují jako neplatné směry způsobené
šumem nebo chybným stavem.

Celkový basis index používá Qiskit little-endian konvenci:

\[
I=x+2^{n_x}y+2^{n_x+n_y}d.
\]

Implementace:

```python
EncodeIndex(x, y, direction, reg_x_size, reg_y_size)
DecodeIndex(index, reg_x_size, reg_y_size)
PrepareBasisState(x, y, direction)
```

Například mřížka 512×512 vyžaduje

```text
9 qubitů x + 9 qubitů y + 4 qubity d = 22 logických qubitů.
```

`ISA qubits` ve výpisu může odpovídat velikosti celého backendu, například
156, i když obvod měří pouze své logické registry.

## Měření a převod counts na pole

Sampler provede pro každý snapshot `SHOTS` měření. Pro bitstring
\((x,y,d)\) s počtem výskytů \(n_{xyd}\) se odhadne pravděpodobnost

\[
\hat P(x,y,d)=\frac{n_{xyd}}{N_{\mathrm{shots}}}.
\]

Pro platné směry se ukládá

\[
f_d(y,x)=\hat P(x,y,d),\qquad d=0,\ldots,8.
\]

Neplatné směry 9–15 nejsou vloženy do `f`, ale jejich celková pravděpodobnost
je uvedena v JSON metadatech. Kvůli tomu může platit

\[
\sum_{d,y,x}f_d(y,x)<1
\]

na hlučném QPU.

Současné `f` je normalizovaná pravděpodobnost měření jednoho zakódovaného
stavu. Není to ještě fyzikální pole hmotnostních populací pro celý CFD
výpočet. Vypočtené `density`, `velocity` a `pressure` jsou proto vhodné pro
vizualizaci a kontrolu experimentu, ne jako hotové aerodynamické řešení.

## Metriky kvality QPU

Každý JSON snapshot obsahuje:

- `expected_state_probability` – pravděpodobnost přesně správného
  \((x,y,d)\);
- `correct_direction_probability` – pravděpodobnost očekávaného směru;
- `valid_direction_probability` – součet všech platných směrů 0–8;
- `invalid_direction_probability` – součet směrů 9–15;
- `wall_time_seconds` – reálný čas přípravy, fronty a měření snapshotu.

Pro správnost streamingu je hlavní první metrika. Vysoké
`valid_direction_probability` samo o sobě neznamená správnou pozici ani
správný směr.

## Časová řada

Parametry:

```python
TIME_STEP = 0.1
SNAPSHOT_COUNT = 40
SHOTS = 1000
```

Fyzikální/lattice čas snapshotu \(k\) je

\[
t_k=k\Delta t,\qquad k=1,\ldots,40,
\]

tedy \(t=0.1,\ldots,4.0\).

Pozor na terminologii:

- **snapshot** je jeden časový snímek;
- **measurement shot** je jedno opakované měření stejného kvantového obvodu;
- jeden QPU snapshot nyní znamená jeden IBM Runtime job s `SHOTS` měřeními.

Při této basis-state demonstraci se každý další snapshot připraví z ideálně
očekávané předchozí pozice. Šum QPU se tedy mezi časovými kroky nekumuluje a
naměřený stav není zpětně vložen jako nový fyzikální stav. Pro plnou simulaci
je nutné propagovat celé pole \(f_i\) přes collision, hranice a streaming.

## Referenční inlet CFD výpočet

`src/lbm_solver.py` implementuje klasický BGK výpočet celého pole
\(f_i(y,x)\). Je oddělený od velikosti QPU testu pomocí
`CFD_CELLS_SIZE`.

Výchozí fyzikální nastavení:

```python
CFD_CELLS_SIZE = 128
PHYSICAL_CELL_SIZE_M = 0.25
PHYSICAL_INLET_VELOCITY_M_S = 1.0
LBM_INLET_VELOCITY = 0.05
LBM_RELAXATION_TIME = 0.6
REFERENCE_DENSITY_KG_M3 = 1.225
REFERENCE_PRESSURE_PA = 101325.0
```

Převod rychlosti je

\[
S_u=\frac{u_{\mathrm{physical}}}{u_{\mathrm{lattice}}}
=\frac{1}{0.05}=20\ \mathrm{m\,s^{-1}}.
\]

Pro velikost buňky \(\Delta x=0.25\ \mathrm{m}\) vychází stabilní interní
časový krok

\[
\Delta t_{\mathrm{LBM}}=\frac{\Delta x}{S_u}
=\frac{0.25}{20}=0.0125\ \mathrm{s}.
\]

Jeden výstupní interval 0.1 s proto obsahuje osm interních LBM kroků.
Nastavit přímo \(\Delta t_{\mathrm{LBM}}=0.1\ \mathrm{s}\) by při rychlosti
1 m/s dalo lattice rychlost 0.4, která je pro weakly-compressible D2Q9 příliš
vysoká.

Použité hranice:

- vlevo equilibrium velocity inlet \(u_x=1\ \mathrm{m/s}\);
- vpravo zero-gradient outlet;
- nahoře a dole periodická hranice.

SI převody jsou

\[
\rho_{\mathrm{SI}}=\rho_{\mathrm{ref}}\rho_{\mathrm{LBM}},
\]

\[
\mathbf{u}_{\mathrm{SI}}=S_u\mathbf{u}_{\mathrm{LBM}},
\]

\[
p_{\mathrm{gauge}}
=\rho_{\mathrm{ref}}S_u^2
\left(p_{\mathrm{LBM}}-\frac13\right),
\qquad
p_{\mathrm{absolute}}=p_{\mathrm{ref}}+p_{\mathrm{gauge}}.
\]

Spuštění:

```sh
make run-cfd
```

Výstupy jsou v `output/cfd_series/`. Každé NPZ obsahuje úplné `f`, SI i
lattice varianty density/velocity/pressure a gauge pressure. Na konci se
automaticky vytvoří GIF a MP4.

## Plný hybridní QLBM prototyp

Hybridní solver je implementován v `src/hybrid_qlbm.py`. Zachovává původní
obecnou bránu `StreamG`, ale na rozdíl od basis-state testu pracuje se všemi
buňkami a všemi devíti D2Q9 populacemi současně.

### Rozdělení klasické a kvantové části

| Operace | Aer hybrid | QPU hybrid |
|---|---|---|
| Inicializace \(f_i\) | klasická | klasická |
| Výpočet \(\rho,\mathbf{u},f_i^{eq}\) | klasický | klasický |
| BGK collision | klasická | klasická |
| Příprava amplitudového stavu | kvantový obvod simulovaný na Aer | kvantový obvod na IBM QPU |
| D2Q9 streaming `StreamG` | kvantový obvod simulovaný na Aer | kvantový obvod na IBM QPU |
| Rekonstrukce \(f_i\) | ze statevector pravděpodobností | z naměřených counts |
| Inlet/outlet conditions | klasické | klasické |
| Density, velocity, pressure | klasické | klasické |

Jeden hybridní časový krok má tento datový tok:

```text
f(t)
  -> klasický BGK collision
  -> f*
  -> amplitudová příprava |psi_f*>
  -> kvantový StreamG
  -> pravděpodobnosti/counts
  -> rekonstrukce f(t + dt)
  -> inlet/outlet boundary conditions
  -> density, velocity, pressure
```

BGK collision před streamingem používá

\[
f_i^*=f_i-\frac{1}{\tau}(f_i-f_i^{eq}).
\]

### Amplitudové kódování úplného pole

Použitý stav:

\[
|\psi_f\rangle
=\sum_{x,y,d}
\sqrt{\frac{f_d(y,x)}{\sum_{q,a,b}f_q(b,a)}}\,|x,y,d\rangle.
\]

Po streamingu dávají ideální pravděpodobnosti:

\[
P(x,y,d)=
\frac{f_d^{\mathrm{streamed}}(y,x)}
{\sum_{q,a,b}f_q(b,a)}.
\]

Normalizační konstanta

\[
M=\sum_{d,y,x}f_d(y,x)
\]

se uloží před spuštěním obvodu. Po získání pravděpodobností se populace
rekonstruují vztahem \(f_d=M P(x,y,d)\). Funkce
`probabilities_to_populations` přijímá pouze fyzikální směry 0–8 a reportuje
jejich celkovou pravděpodobnost. Na ideálním Aer běhu musí být tato hodnota
prakticky 100 %. Na QPU její pokles ukazuje únik do neplatných direction
stavů 9–15.

### Kvantový streaming

`StreamG` realizuje unitární periodický posun

\[
|x,y,d\rangle\mapsto
|(x+c_{d,x})\bmod N_x,\,
  (y+c_{d,y})\bmod N_y,\,
  d\rangle.
\]

Protože je operace lineární, jedna aplikace brány streamuje amplitudy všech
buněk a směrů v superpozici. Test
`tests/test_hybrid_qlbm.py` porovnává rekonstruované celé pole s klasickým
D2Q9 streamingem, nikoliv pouze jeden basis state.

### Aer režim

Aer varianta používá `AerSimulator(method="statevector")`. Pravděpodobnosti
se získají z přesného statevectoru, takže zde nevzniká sampling noise:

```sh
make run-hybrid-aer
```

Výchozí konfigurace provede 40 časových kroků po 0,1 s, tedy interval
0–4 s. Výstupy jsou v:

```text
output/hybrid_aer/
```

### Experimentální QPU režim

Stejnou amplitudovou přípravu a `StreamG` lze transpileovat a spustit na IBM
QPU:

```sh
make run-hybrid-qpu
```

Výstupy:

```text
output/hybrid_qpu/
```

QPU vrací konečný počet měření. Pro bitstring odpovídající \((x,y,d)\)
odhadujeme

\[
\hat P(x,y,d)=\frac{n_{xyd}}{N_{\mathrm{shots}}}.
\]

Rekonstrukce úplného pole je proto zatížena sampling noise, gate errors,
readout errors a chybami z hlubokého transpileovaného obvodu. Každý další
časový krok navíc používá rekonstruované hlučné pole jako vstup následujícího
kroku, takže chyby se mohou kumulovat.

Výchozí hybridní síť je z tohoto důvodu pouze 4×4. Aer provádí 40 kroků,
zatímco QPU režim standardně provede jediný experimentální krok:

```python
HYBRID_CELLS_SIZE = 4
HYBRID_LATTICE_INLET_VELOCITY = 0.03
HYBRID_RELAXATION_TIME = 0.8
HYBRID_SNAPSHOT_COUNT = 40
HYBRID_QPU_SNAPSHOT_COUNT = 20
HYBRID_SHOTS = 32768
```

Při `TIME_STEP = 0.1` odpovídá 20 QPU snímků časovému intervalu
\(0.1,\ldots,2.0\) s. Každý časový krok vytváří samostatný QPU job se
`HYBRID_SHOTS = 32768` měřeními; celá série tedy vyžaduje 20 hardwarových
jobů a až 655 360 měření. Pole rekonstruované z jednoho hlučného kroku se
použije jako vstup dalšího kroku, takže se hardwarové chyby v čase kumulují.

Výchozí QPU režim není validované aerodynamické řešení. Obecná příprava
amplitud a direction-controlled streaming vytvářejí hluboký obvod. Při
hodnocení experimentu je proto nutné kontrolovat:

- `valid_direction_probability`;
- hloubku transpileovaného obvodu;
- počet dvouqubitových `cz` bran;
- zachování celkové hmotnosti;
- rozsah density, velocity a pressure;
- odchylku od ideálního Aer výsledku.

### Výstupní data

Každý krok ukládá numerická data a vizualizaci. Výstupní adresář obsahuje
zejména:

```text
step_XXXX.npz          úplné f a makroskopická pole
step_XXXX.json         metadata a diagnostika běhu
step_XXXX.png          density, pressure a velocity
manifest.json          seznam a parametry kroků
qlbm_animation.gif     časová animace
qlbm_animation.mp4     časová animace
```

NPZ soubory lze načíst například takto:

```python
import numpy as np

data = np.load("output/hybrid_aer/step_0040.npz")
f = data["f"]
density = data["density"]
velocity = data["velocity"]
pressure = data["pressure"]
```

### Co hybridní QLBM znamená

Kvantová část nyní nahrazuje streamingový operátor klasického LBM. Neznamená
to automatické kvantové zrychlení celého CFD výpočtu: příprava obecného
amplitudového stavu a získání celého pole měřením jsou nákladné a QPU hardware
je hlučný. Implementace slouží jako korektní funkční prototyp a základ pro
další práci na efektivnější přípravě stavů, collision obvodu, error mitigation
a fyzikálních boundary conditions včetně překážky.

## Spuštění

Instalace:

```sh
make i
```

Jeden lokální snapshot:

```sh
make run-local
```

Jeden QPU snapshot:

```sh
make run-qpu
```

Časová série na lokálním simulátoru:

```sh
make run-series-local
```

Časová série na IBM QPU:

```sh
make run-series-qpu
```

Plný referenční BGK inlet výpočet:

```sh
make run-cfd
```

Plný hybridní 4×4 výpočet na Aer:

```sh
make run-hybrid-aer
```

Jeden experimentální plný hybridní krok na QPU:

```sh
make run-hybrid-qpu
```

Oba time-series cíle nejprve odstraní `output/time_series/`, aby se nemíchala
stará a nová měření. QPU varianta odešle `SNAPSHOT_COUNT` samostatných Runtime
jobů.

Testy:

```sh
make test-qlbm
```

Zobrazení kvantových obvodů:

```sh
make show-circuits
```

## Výstupní data a animace

Každý snímek v `output/time_series/` vytvoří:

- `step_NNNN.png` – vizualizaci prostorové distribuce, směrů a rychlosti;
- `step_NNNN.npz` – numerická pole;
- `step_NNNN.json` – metadata a metriky QPU;
- `step_NNNN.csv` – surové dekódované counts.

Společné soubory:

- `manifest.json` – časování a seznam všech snapshotů;
- `qlbm_animation.gif` – animace s periodou `TIME_STEP`;
- `qlbm_animation.mp4` – kompaktní video, pokud je dostupný `ffmpeg`.

Animaci lze znovu vytvořit bez nového měření:

```sh
make animation
```

Načtení numerických dat:

```python
import numpy as np

snapshot = np.load("output/time_series/step_0001.npz")
f = snapshot["f"]                # (9, ny, nx)
density = snapshot["density"]    # (ny, nx)
velocity = snapshot["velocity"]  # (ny, nx, 2)
pressure = snapshot["pressure"]  # (ny, nx)
```

## Okrajové podmínky a další vývoj

Hybridní prototyp obsahuje úplné pole populací, BGK collision, inlet/outlet
a převod makroskopických polí do SI jednotek. Samotný `StreamG` používá
periodické modulo; po rekonstrukci pole se jeho periodické chování v ose x
přepíše dvěma explicitními hranicemi:

- `x = 0`: pevný rychlostní inlet s \(\rho=1\),
  \(u_x=u_{\mathrm{inlet}}\) a \(u_y=0\);
- `x = N_x-1`: pravý zero-gradient outlet, kde
  \(f_i(x=N_x-1)=f_i(x=N_x-2)\).

Horní a dolní okraj zůstávají periodické. Stejnou funkci
`apply_channel_boundaries` používá klasický i hybridní solver, takže jejich
okrajové podmínky jsou přímo porovnatelné.

Pevný inletový sloupec `x=0` se ukládá do NPZ souborů, ale nezobrazuje se
v PNG grafech, v animacích ani se nezapočítává do časových diagnostik.
Vizualizace proto začínají první vnitřní buňkou `x=1`; numerický stav solveru
zůstává kompletní a nezměněný.

Pro validovanou simulaci obtékání křídla je ještě potřeba doplnit:

1. bounce-back nebo interpolovanou no-slip podmínku na povrchu křídla;
2. masku pevné oblasti `solid_mask.shape == (ny, nx)`;
3. konzistentní volbu rozlišení, \(\tau\), viskozity a Reynoldsova čísla;
4. aerodynamické síly, koeficienty vztlaku a odporu;
5. validaci proti klasickému LBM nebo OpenFOAM;
6. QPU error mitigation a snížení hloubky obvodu.

Dokud tyto části a validační testy nejsou dokončené, nelze QPU tlak
interpretovat jako spolehlivý tlak na křídle ani z něj určovat vztlak a odpor.

## IBM credentials

Lokální soubor `src/credentials.py`:

```python
TOKEN = "..."
CRN = "..."
```

Credentials necommitujte do Git repozitáře.
