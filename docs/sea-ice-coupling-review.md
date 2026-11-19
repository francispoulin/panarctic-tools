# Sea-ice coupling: our model against GLORYS12 (NEMO 3.1 + LIM2), and how NEMO's sea ice is tested

19 November 2026. Context: PanArcticShenanigans issue #48 (surface layer too salty under ice; CATKE and RiBased alike,
much better with no closure).

Code read: PanArcticShenanigans `main` (fc38201); NumericalEarth `fjp/net-inflow` (the pinned #938 on #676);
ClimaSeaIce v0.7.2 (the Manifest's version). NEMO side: GLORYS12's published configuration (QUID, technical sheet,
Lellouche et al. 2021), plus LIM2 and TKE source from a NEMO 3.6-era tree (GLORYS12 runs NEMO 3.1; the coupling
formulas are believed to be the same, but GLORYS's namelist values are not published).

## 1. Main findings

1. **The ice–ocean stress is coupled the same way in both models**: quadratic drag (ours 5.5e-3, ORCA-LIM 5e-3),
   against the top ocean cell's velocity, no turning angle (LIM2 forces it to zero with EVP), stress on the ocean
   weighted by ice concentration. Our top cell is 8 m, GLORYS's 1 m. So the coupling formula is not the difference.
2. **Neither GLORYS's NEMO nor CATKE reduces the surface TKE input under ice, but NEMO damps everything else.** In
   NEMO's TKE scheme the surface value is `rn_ebb × τ/ρ`, with τ including the ice–ocean stress, as in CATKE. But the
   extra mixing NEMO adds below the mixed layer (`nn_etau`, "ETAU", e-folding scale 10 m in GLORYS12) and the
   Langmuir source are multiplied by (1 − ice fraction). NEMO 4 added `nn_eice`, which also attenuates the surface
   wave-breaking input under ice (default 1 − tanh(10 fᵢ)).
3. **The literature has our symptom.** Allende et al. (2024, GMD) ran NEMO-SI3 in the Arctic basins: without the
   ice scaling of TKE penetration, the Canada Basin surface became more than 2 psu saltier, mixed layers deepened by
   20 m or more, and stratification weakened. Their conclusion: the TKE input must be scaled down under sea ice.
   CATKE has no such scaling. This is the strongest lead so far, and it is "using CATKE properly under ice" rather
   than tuning it.
4. **Ice–ocean heat exchange is 2–6× weaker in ours.** Ours: three-equation formulation with a **constant** u★ =
   0.002 m/s and Γ_T = 0.0095 (exchange velocity 1.9e-5 m/s), whatever the stress. LIM2: McPhee bulk, 0.006 × u★, with
   u★ from the actual ice–ocean stress clamped to 0.5–2 cm/s (3e-5 to 1.2e-4 m/s). Weaker exchange keeps heat in the
   ocean under ice and slows basal melt (less fresh water at the surface).
5. **Leads close 20× faster in ours.** Open-water freezing raises concentration as (1 − ℵ)/h₀ with h₀ = the
   consolidation thickness, 0.05 m (ClimaSeaIce default via NumericalEarth); LIM2 uses h₀ = 1.0 m (ORCA-LIM), Hibler
   0.5 m. With 0.05 m, freeze-up turns open water into near-full ice cover almost at once: less open-water heat loss,
   less frazil spread, and the ocean's stress switches from wind to ice stress quickly. Matters most from October.
6. **GLORYS is not a reliable target at the Arctic surface.** Its own QUID reports a salty bias in the top 500 m of
   the Beaufort Sea and Arctic surface values saltier than climatology (wording inconsistent within the document);
   its ice concentration is assimilated (with thickness adjusted), and assimilation is said to weaken leads and
   thicken ice. Agreement with GLORYS at the surface is a consistency check, not validation.

## 2. Side by side

| Component | Ours | GLORYS12 (NEMO 3.1 + LIM2) | Possible consequence |
|---|---|---|---|
| Ice thermodynamics | Zero-layer slab (no heat capacity), slab snow; ice salinity 4 | 2 ice layers + 1 snow layer with heat storage, brine reservoir; ice salinity 6 (ORCA-LIM) | Faster melt/freeze response; ~7 % more salt rejected per metre frozen |
| Categories | One (h, ℵ) | One (LIM2) | — |
| Rheology | EVP (adaptive, 150 substeps), P* = 27 500 N/m², C = 20, e = 2 | EVP; P* 20 000 in ORCA-LIM (GLORYS value not published) | Similar |
| Ice–ocean drag | 5.5e-3, top cell (8 m), no turning, semi-implicit both sides | 5e-3, top level (1 m), no turning | Same formula; different cell depth |
| Ice–ocean heat flux | Three-equation, **constant u★ 0.002**, Γ_T 0.0095 | McPhee bulk 0.006 u★, u★ from stress, clamped 0.5–2 cm/s | Ours 2–6× weaker ocean→ice heat flux |
| Salt/fresh water | Ice melt/growth as a volume flux; virtual salt flux with local surface S; ice salt Sᵢ | Virtual salt flux, levitating ice (same idea) | Same |
| Lead closing h₀ | 0.05 m (= consolidation thickness) | 1.0 m (ORCA-LIM) | Ours closes leads ~20× faster at freeze-up |
| Albedo | Constant 0.7 (`SurfaceRadiationProperties(0.7, 1.0)` default) | Shine & Henderson-Sellers: depends on melting, snow, thickness | Ours too bright for melting bare ice (~0.5–0.6) |
| Air–ice fluxes | Similarity theory, ocean (Charnock) roughness with ice stability functions | CORE (Large & Yeager) constant Cᵢ = 1.4e-3 | Neutral Cd ≈ 1.35e-3 in ours at typical winds: similar or lower wind stress |
| Atmosphere | JRA55 | ERA-Interim, radiation not corrected north of 65°N | Different forcing (ERA5 PR planned) |
| Vertical mixing | CATKE; surface TKE from total stress incl. ice | TKE 1.5; surface from total stress incl. ice; ETAU and Langmuir × (1 − fᵢ); background 1e-5 / 1e-4 | No under-ice damping in ours |
| Top cell | ~8 m (nz50), ~4 m (nz100) | ~1 m (22 levels in top 100 m) | Ice drag and brine into a thick cell |
| Initial ice | GLORYS h, ℵ; **snow 0, ice velocity 0** | Spun up, assimilated | Spin-up of drift and snow insulation |
| Rivers | JRA55 runoff | Dai et al. (2009) climatology | Different shelf freshwater |

## 3. What they have that we don't

- Under-ice damping of TKE input (ETAU and Langmuir in NEMO 3; also surface wave breaking via `nn_eice` in NEMO 4).
- Stress-dependent ice–ocean heat exchange (our u★ is constant unless `friction_velocity =
  MomentumBasedFrictionVelocity()` is passed to `ThreeEquationHeatFlux`, which exists in NumericalEarth).
- A lead-closing parameter separate from the consolidation thickness.
- State-dependent albedo (NumericalEarth has `SeaIceAlbedo` (CCSM3, Briegleb et al. 2004) but it is not the default).
- Thermodynamic heat capacity (layers, brine). Not cheap to add; lower priority.
- Assimilation of concentration (not wanted for a free model, but it explains part of the GLORYS difference).

## 4. How NEMO's sea ice (and coupled sea ice in general) is tested

**NEMO test cases** (`tests/` in NEMO 4): `ICE_ADV1D` (Schär convergence test: consistency of concentration,
thickness and volume), `ICE_ADV2D` (advection schemes and overshoots), `ICE_AGRIF` (ice across a nest boundary),
`ICE_RHEO` (rheology), `C1D_ASICS` (a 1-D column configuration), `STATION_ASF` (bulk air–sea flux formulae at a
station, pass/fail). These are component tests; none validates the coupled Arctic.

**Coupled-model evaluation** uses observations, not test cases:
- SITool (Lin et al. 2021, GMD; used for CMIP6 OMIP): concentration and extent (NSIDC-0051, OSI-450), ice edge
  (integrated ice-edge error), thickness (Envisat, ICESat), snow depth (Envisat, SnowModel-LG), drift magnitude and
  direction (NSIDC v4.1 ×1.357, KIMURA). Each error is scaled by the difference between two observational products,
  so 1 means "as good as the observations agree".
- GLORYS12's own validation: extent and concentration (CERSAT, NSIDC, OSI SAF), thickness (CryoSat-2, BGEP moorings,
  NPEO), volume (PIOMAS); ocean T/S against CORA profiles.
- ORCA-LIM (Timmermann et al. 2005): drift (SSM/I, QuikSCAT), Fram Strait ice export, submarine thickness, PHC
  hydrography.
- Arctic upper ocean: Ice-Tethered Profilers (as in Allende et al. 2024), BGEP moorings, WOA/PHC climatologies.

## 5. Proposed tests

In order, each one change against the 30-day CATKE baseline, using `salinity_drift.jl` (Beaufort, Central Arctic):

1. **Under-ice scaling of CATKE's surface TKE flux**, (1 − ℵ) and 1 − tanh(10 ℵ). The scale-0 run (patch 070) is
   the upper bound.
2. **Stress-dependent ice–ocean heat flux**: `ThreeEquationHeatFlux(sea_ice; friction_velocity =
   MomentumBasedFrictionVelocity())`. Watch basal melt and surface S; check u★ stays reasonable (fast ice).
3. **Lead closing**: a separate h₀ (0.5–1 m) for the concentration step. Probably needs a ClimaSeaIce option; most
   relevant for October onward.
4. **Albedo**: `SeaIceAlbedo` instead of 0.7. Small in September, larger for a melt season.
5. **Initial ice velocity and snow** from GLORYS if available; otherwise note the spin-up.

And for validation beyond GLORYS: SITool-style concentration/extent/drift metrics (NSIDC, OSI SAF), ITP profiles in
the Beaufort for Sep–Oct 2005. Patch 069 (ice speed vs GLORYS) is the first step for drift.

## Sources

- Lellouche et al. (2021), The Copernicus Global 1/12° Oceanic and Sea Ice GLORYS12 Reanalysis, Front. Earth Sci.
  https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2021.698876/full
- CMEMS GLO QUID 001-030: https://documentation.marine.copernicus.eu/QUID/CMEMS-GLO-QUID-001-030.pdf
- Mercator GLORYS12 technical sheet: https://www.mercator-ocean.eu/wp-content/uploads/products_pdf/MOI_GLORYS12.pdf
- Timmermann et al. (2005), ORCA-LIM: https://epic2-clone.awi.de/10070/1/Tim2004b.pdf
- NEMO 3.6-era LIM2/TKE source (limsbc_2, limthd_2, zdftke, sbcblk_core): https://github.com/fmassonn/paper-itd-seaice
- NEMO `namzdf_tke` (nn_eice): https://svn-mirror.nemo-ocean.eu/NEMO/trunk/doc/namelists/namzdf_tke
- Allende et al. (2024), GMD 17, 7445: https://gmd.copernicus.org/articles/17/7445/2024/
- NEMO test cases: https://forge.nemo-ocean.eu/hatfield/nemo/-/raw/e6c54e2b21ba103f08095d4395e3a728f750e224/tests/README.rst
- Lin et al. (2021), SITool, GMD 14, 6331: https://gmd.copernicus.org/articles/14/6331/2021/
