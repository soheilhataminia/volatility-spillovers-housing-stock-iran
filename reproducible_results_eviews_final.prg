'============================================================
' Reproducible EViews Program: Housing–Stock Volatility (Iran)
' Period: 2016M04–2023M10
' Pipeline: ARMA -> GARCH(1,1) -> DCC(1,1) + Cheung–Ng (m=12) + PP
' EViews: 13 or later (MGARCH DCC supported)
'============================================================

'------------------------------
' 0) Workfile & Data Import
'------------------------------
wfcreate(wf=MGARCH_IR, page=Main) m 2016m04 2023m10

' Excel file must sit next to this PRG:
%datapath = @runpath + "\mgarch_gregorian_unmodified.xlsx"

' Import from sheet "gregorian" with columns:
' date_gregorian (YYYY-MM), stock_index_close, tehran_price_per_m2_rial
import(options) {%datapath} range=gregorian @date(date_gregorian,"YYYY-MM") stock_index_close tehran_price_per_m2_rial

rename stock_index_close capidx
rename tehran_price_per_m2_rial estp

'------------------------------
' 1) Returns (log * 100)
'------------------------------
series rcap = 100*(log(capidx) - log(capidx(-1)))
series rest = 100*(log(estp)  - log(estp(-1)))

' Use only non-missing observations
smpl @all
smpl if @isna(rcap)=0 and @isna(rest)=0

'------------------------------
' 2) Descriptives, ADF, PP
'------------------------------
spool sp_desc
spool sp_adf
spool sp_pp

sp_desc.append rcap.stats
sp_desc.append rest.stats

equation adf_rcap.adf rcap 0 0 c
equation adf_rest.adf rest 0 0 c
sp_adf.append adf_rcap.output
sp_adf.append adf_rest.output

' Phillips–Perron (level, intercept only) for Appendix A5
freeze(pp_rcap)  rcap.unitroot(pp, c)
freeze(pp_rest)  rest.unitroot(pp, c)
sp_pp.append pp_rcap
sp_pp.append pp_rest

'------------------------------
' 3) ARMA order selection (AIC; p,q in {0,1,2})
'------------------------------
subroutine buildspec(string %y, scalar !p, scalar !q) string %spec
  %spec = %y + " c"
  if !p>0 then
    for !i=1 to !p
      %spec = %spec + " ar(" + @str(!i) + ")"
    next
  endif
  if !q>0 then
    for !j=1 to !q
      %spec = %spec + " ma(" + @str(!j) + ")"
    next
  endif
endsub

' RCAP
scalar !best_aic_rc = 1e+30
scalar !best_p_rc  = 0
scalar !best_q_rc  = 0
for !p=0 to 2
  for !q=0 to 2
    call buildspec("rcap", !p, !q) %s
    equation eqtmp1.ls {%s}
    if @aic < !best_aic_rc then
      !best_aic_rc = @aic
      !best_p_rc = !p
      !best_q_rc = !q
      eq_rcap_mean.copy eqtmp1
    endif
  next
next

' REST
scalar !best_aic_re = 1e+30
scalar !best_p_re  = 0
scalar !best_q_re  = 0
for !p=0 to 2
  for !q=0 to 2
    call buildspec("rest", !p, !q) %s
    equation eqtmp2.ls {%s}
    if @aic < !best_aic_re then
      !best_aic_re = @aic
      !best_p_re = !p
      !best_q_re = !q
      eq_rest_mean.copy eqtmp2
    endif
  next
next

' ARMA order table
table tbl_arma
tbl_arma.setwidth(1) 24
tbl_arma.setwidth(2) 22
tbl_arma.setwidth(3) 14
tbl_arma.setwidth(4) 14
tbl_arma.setelem(1,1) Series
tbl_arma.setelem(1,2) Selected ARMA(p,q)
tbl_arma.setelem(1,3) AIC
tbl_arma.setelem(1,4) BIC
scalar __bic_rc = eq_rcap_mean.@schwarz
scalar __bic_re = eq_rest_mean.@schwarz
tbl_arma.setelem(2,1) RCAP
tbl_arma.setelem(2,2) ({!best_p_rc}, {!best_q_rc})
tbl_arma.setelem(2,3) {!best_aic_rc}
tbl_arma.setelem(2,4) {__bic_rc}
tbl_arma.setelem(3,1) REST
tbl_arma.setelem(3,2) ({!best_p_re}, {!best_q_re})
tbl_arma.setelem(3,3) {!best_aic_re}
tbl_arma.setelem(3,4) {__bic_re}

'------------------------------
' 4) Residual diagnostics for ARMA (12 lags)
'------------------------------
series e_stock = eq_rcap_mean.resid
series e_house = eq_rest_mean.resid

freeze(lb_resid_rc)  e_stock.qstat(12)
freeze(lb_resid_re)  e_house.qstat(12)
freeze(lb_sq_rc)     (e_stock^2).qstat(12)
freeze(lb_sq_re)     (e_house^2).qstat(12)

freeze(archlm_rc)    eq_rcap_mean.archlm(12)
freeze(archlm_re)    eq_rest_mean.archlm(12)

'------------------------------
' 5) Univariate GARCH(1,1) on ARMA residuals (mean=0)
'------------------------------
equation g_stock.arch e_stock c arch(1) garch(1)
equation g_house.arch e_house c arch(1) garch(1)

series h_stock = g_stock.@h
series h_house = g_house.@h
series s_stock = @sqrt(h_stock)
series s_house = @sqrt(h_house)
series z_stock = e_stock / s_stock
series z_house = e_house / s_house

'------------------------------
' 6) DCC(1,1) estimation and dynamic correlation
'------------------------------
group GZ z_stock z_house
equation dcc_eq.mgarch(dcc, garch=(1,1), dist=normal, type=diagonal) z_stock z_house
freeze(dcc_out) dcc_eq.output

' Reconstruct DCC-implied correlation series (rho_t)
scalar a = @coef("dcc_eq","dcc_a")
scalar b = @coef("dcc_eq","dcc_b")
scalar s11 = @var(z_stock)
scalar s22 = @var(z_house)
scalar s12 = @cov(z_stock,z_house)

series q11 = s11
series q22 = s22
series q12 = s12
smpl @first+1 @last
q11 = (1-a-b)*s11 + a*z_stock(-1)^2            + b*q11(-1)
q22 = (1-a-b)*s22 + a*z_house(-1)^2            + b*q22(-1)
q12 = (1-a-b)*s12 + a*z_stock(-1)*z_house(-1)  + b*q12(-1)
series rho_t = q12/@sqrt(q11*q22)
smpl @all

' Summary stats for rho_t
scalar rho_mean = @mean(rho_t)
scalar rho_min  = @min(rho_t)
scalar rho_max  = @max(rho_t)
scalar rho_sd   = @stdev(rho_t)

table tbl_dcc
tbl_dcc.setwidth(1) 20
tbl_dcc.setwidth(2) 16
tbl_dcc.setelem(1,1) Statistic
tbl_dcc.setelem(1,2) Value
tbl_dcc.setelem(2,1) Mean
tbl_dcc.setelem(2,2) {rho_mean}
tbl_dcc.setelem(3,1) Min
tbl_dcc.setelem(3,2) {rho_min}
tbl_dcc.setelem(4,1) Max
tbl_dcc.setelem(4,2) {rho_max}
tbl_dcc.setelem(5,1) Std. Dev.
tbl_dcc.setelem(5,2) {rho_sd}

'------------------------------
' 7) Cheung–Ng (1996) variance-causality, m=12
'------------------------------
series u = z_stock^2 - @mean(z_stock^2)
series v = z_house^2 - @mean(z_house^2)

vector(12) rho_sh
vector(12) rho_hs
scalar Tobs = @obs(u)

for !k=1 to 12
  smpl @first+!k @last
  series u_lag = u(-!k)
  series v_lag = v(-!k)
  rho_sh(!k) = @cor(v, u_lag)   ' Stock -> Housing
  rho_hs(!k) = @cor(u, v_lag)   ' Housing -> Stock
next
smpl @all

scalar Q_sh = Tobs * @inner(rho_sh, rho_sh)
scalar Q_hs = Tobs * @inner(rho_hs, rho_hs)
scalar p_sh = 1 - @cchisq(Q_sh, 12)
scalar p_hs = 1 - @cchisq(Q_hs, 12)

table tbl_cn
tbl_cn.setwidth(1) 24
tbl_cn.setwidth(2) 20
tbl_cn.setwidth(3) 16
tbl_cn.setelem(1,1) Direction
tbl_cn.setelem(1,2) Q-stat (df=12)
tbl_cn.setelem(1,3) p-value
tbl_cn.setelem(2,1) Stock -> Housing
tbl_cn.setelem(2,2) {Q_sh}
tbl_cn.setelem(2,3) {p_sh}
tbl_cn.setelem(3,1) Housing -> Stock
tbl_cn.setelem(3,2) {Q_hs}
tbl_cn.setelem(3,3) {p_hs}

table tbl_cn_lag
tbl_cn_lag.setwidth(1) 6
tbl_cn_lag.setwidth(2) 28
tbl_cn_lag.setwidth(3) 28
tbl_cn_lag.setelem(1,1) Lag
tbl_cn_lag.setelem(1,2) Corr[stock^2→housing]
tbl_cn_lag.setelem(1,3) Corr[housing^2→stock]
for !k=1 to 12
  tbl_cn_lag.setelem(1+!k,1) {!k}
  tbl_cn_lag.setelem(1+!k,2) {rho_sh(!k)}
  tbl_cn_lag.setelem(1+!k,3) {rho_hs(!k)}
next

'------------------------------
' 8) Figures: Conditional Volatilities & DCC Correlation
'   (Renumbered to match paper: these are Figures 2 and 3)
'------------------------------
graph fig_vol.line s_stock s_house
fig_vol.addtext(t) "Figure 2. Conditional Volatilities (GARCH(1,1))"
fig_vol.legend on
fig_vol.axis(l) label "Sigma"
fig_vol.save(t=png) @runpath + "\Figure2_conditional_vols.png"

graph fig_dcc.line rho_t
fig_dcc.addtext(t) "Figure 3. DCC Dynamic Correlation: Stock vs Housing"
fig_dcc.legend off
fig_dcc.axis(l) label "Correlation"
fig_dcc.save(t=png) @runpath + "\Figure3_dcc_correlation.png"

'------------------------------
' 9) Spool & Export
'------------------------------
spool SP_ALL
SP_ALL.append sp_desc
SP_ALL.append sp_adf
SP_ALL.append sp_pp
SP_ALL.append eq_rcap_mean.output
SP_ALL.append eq_rest_mean.output
SP_ALL.append tbl_arma
SP_ALL.append lb_resid_rc
SP_ALL.append lb_resid_re
SP_ALL.append lb_sq_rc
SP_ALL.append lb_sq_re
SP_ALL.append archlm_rc
SP_ALL.append archlm_re
SP_ALL.append g_stock.output
SP_ALL.append g_house.output
SP_ALL.append dcc_out
SP_ALL.append tbl_dcc
SP_ALL.append tbl_cn
SP_ALL.append tbl_cn_lag
SP_ALL.append fig_vol
SP_ALL.append fig_dcc

' Export tables/views to Excel
export(table, mode=overwrite) @runpath + "\results_tables_from_eviews.xlsx" tbl_arma tbl_dcc tbl_cn tbl_cn_lag pp_rcap pp_rest

' Save workfile
save @runpath + "\MGARCH_IR_EViews.wf1"

'================== End of Program ==================
