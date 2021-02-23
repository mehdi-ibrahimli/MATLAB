function [H,a,URL] = load_data(fname)

% [H,a,URL] = load_mv(filename, verbose)
%
% loads a .dat file of webgraph
%
% filename is the file name of the .dat file
%
% output 
%  H is the hyperlink matrix where the diagonal Elements are all zero 
%    (no self links). The data format is sparse.
%
%  a is the vector where a_i is 1 if site i contains no link and otherwise 0
%
%  URL is a cell array of the pages urls (unified resource identifier),
%      which are the 'http://www.domain.de/directory/file.html' like strings.
%
% examples
%  [H,a,URL] = load_mv('na_math_kit.dat');
%  [H,a,URL] = load_mv('math_kit.dat');


file   = fopen(fname,'r');
dim    = fscanf(file,'%d',1); % dimension (number of pages)
nnzero = fscanf(file,'%d',1); % number of nonzero entries (links)
fclose(file);

% read data block of urls
[URL]   = textread(fname,'%*d%s',dim,'headerlines',1);
% read data block of web graph (adjacence matrix)
[ii,jj] = textread(fname,'%d%d', nnzero,'headerlines',1+dim);
% size(ii)
% size(jj)
% nnzero

%create adjacence matrix
% put the date into a sparse matrix format
A       = sparse(ii,jj,ones(nnzero,1),dim,dim);
% remove self linking
A       = A-spdiags(diag(A), 0, dim,dim); 

%Create vector a and hyperlink matrix H  
s = sum(A,2);
a = (s==0);
s = s + a; % to avoid division by zero (0/0)
H = spdiags(1./s, 0, dim,dim) * A; % save memory via loop here
clear('A');


return;