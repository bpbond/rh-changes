



srdb_filtered %>% filter(FLUXNET_DIST < 1) -> x

ggplot(x, aes(Year, Rs_annual/GPP_NT_VUT_REF, group=FLUXNET_SITE_ID)) + geom_point(aes(color=Stage)) + facet_grid(sign(tmp_trend) ~ sign(pre_trend), scales="free") + geom_smooth(method='lm', color='grey', fill=NA) + geom_smooth(method="lm", group=1)


