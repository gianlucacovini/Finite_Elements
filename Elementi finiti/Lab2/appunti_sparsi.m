clear
close all

% A = spalloc(n_v, n_v, n_nz)
% n_v = dimensioni, n_nz = numero di non zero
% mettiamo come n_nz almeno 7*n_v, ma meglio essere larghi per non finire
% out of memory, quindi mettiamo 8*n_v
% sparse è più efficiente se poi costruiamo blocco per blocco (il tempo è
% lineare e non superlineare rispetto a n_nz), però ci sono altre
% complicazioni.
