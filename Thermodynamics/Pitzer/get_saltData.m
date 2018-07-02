function W_g_mol = get_saltData(m_molal)

W_g_mol = zeros(1,length(m_molal));
ion_str = fieldnames(m_molal);
n_ions = length(ion_str);
for hi = 1:length(m_molal)
    for hii = 1:n_ions
            m_this = getfield(m_molal(hi),ion_str{hii});
        switch ion_str{hii};
            case 'Na'
            W_g_mol(hi) = W_g_mol(hi)+m_this*22.99;
            case 'Mg'
            W_g_mol(hi) = W_g_mol(hi)+m_this*24.31;
            case 'Ca'
            W_g_mol(hi) = W_g_mol(hi)+m_this*40.08;
            case 'Cl'
            W_g_mol(hi) = W_g_mol(hi)+m_this*35.4527;
            case 'SO4'
            W_g_mol(hi) = W_g_mol(hi)+m_this*(32.06 + 4*16.00);
        end
    end
end