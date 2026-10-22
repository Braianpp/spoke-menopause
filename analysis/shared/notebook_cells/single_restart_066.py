# Credential-free source extraction. See docs/SOURCE_MAP.csv.
import os
personalization = {node_id: 1}
try:
    pagerank = nx.pagerank(G, alpha=0.85, personalization=personalization, max_iter=500, tol=1e-20)
    embedding = np.array(list(pagerank.values()))
    features = np.array(list(pagerank.keys()))
    embedding_df = pd.DataFrame(zip(features, embedding), columns=['node', 'pagerank'])
    top_pagerank_nodes = embedding_df.sort_values(by='pagerank', ascending=False).head(20)
except Exception as e:
    print(e)
