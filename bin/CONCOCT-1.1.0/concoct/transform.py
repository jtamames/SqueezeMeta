from sklearn.decomposition import PCA
import pandas as p

def perform_pca(d, nc, seed):
    pca_object = PCA(n_components=nc, random_state=seed).fit(d)
    return pca_object.transform(d), pca_object

def perform_split_pca(cov_d, composition, pca_components, use_pcas = None):
    if use_pcas:
        cov_pca, comp_pca = use_pcas
        tf_cov = cov_pca.transform(cov_d)
        tf_comp = comp_pca.transform(composition)
    else:
        tf_cov, cov_pca = perform_pca(cov_d, pca_components[0])
        tf_comp, comp_pca = perform_pca(composition, 
                                        pca_components[1])

    tf_cov = p.DataFrame(tf_cov, index=cov_d.index)
    tf_cov = tf_cov.rename(columns=lambda x: 'cov_'+str(x))

    tf_comp = p.DataFrame(tf_comp,
                          index=composition.index)
    tf_comp = tf_comp.rename(columns=lambda x: 'comp_'+str(x))

    return tf_comp.join(tf_cov, how='inner'), cov_pca, comp_pca

