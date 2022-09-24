import pickle

vocab = pickle.load(open('vocab.pkl', 'rb'))

print(vocab.w2i)
