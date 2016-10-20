function I = get_ionic_strength(m_molal)
I = 0;
if isfield(m_molal,'Ca')            
    z_cat_sq = 2;
    I = calc_I(z_cat_sq,[m_molal.Ca],I);
end
if isfield(m_molal,'Na')            
    z_cat_sq = 1;
    I = calc_I(z_cat_sq,[m_molal.Na],I);
end
if isfield(m_molal,'Mg')
    z_cat_sq = 4;
    I = calc_I(z_cat_sq,[m_molal.Mg],I);
end
if isfield(m_molal,'Fe')
    z_cat_sq = 4;
    I = calc_I(z_cat_sq,[m_molal.Fe],I);
end
if isfield(m_molal,'K')
    z_cat_sq = 1;
    I = calc_I(z_cat_sq,[m_molal.K],I);
end
if isfield(m_molal,'Cl')
    z_an_sq = 1;
    I = calc_I(z_an_sq,[m_molal.Cl],I);
end
if isfield(m_molal,'ClO3')
    z_an_sq = 1;
    I = calc_I(z_an_sq,[m_molal.ClO3],I);
end
if isfield(m_molal,'SO4')    
    z_an_sq = 4;
    I = calc_I(z_an_sq,[m_molal.SO4],I);
end

function I = calc_I(z_sq,molal_vec,I)
    I = I+0.5*([molal_vec]*z_sq);
