function createnips
% produces corpus format for knowceans-corpus/tools
% from gal chechnik's nips1-17 data at
% http://ai.stanford.edu/~gal/Data/NIPS/nips_1-17.mat

load 'nips_1-17';
mkdir('nips17');

% term-document matrix
f = fopen('nips17/nips17.corpus', 'w');
for i = 1:size(counts, 2)
    % nwords x ndocs sparse, convert to zero offset
    [terms, freqs] = find(counts(:, i));
    terms = terms - 1;
    fprintf(f, '%d', length(terms));
    for j = 1:length(terms)
        fprintf(f, ' %d:%d', terms(j), freqs(j));
    end
    fprintf(f, '\n');
end
fclose(f);

% vocabulary
f = fopen('nips17/nips17.vocab', 'w');
for i = 1:length(words)
    fprintf(f, '%s\n', words{i});
end
fclose(f);

% document metadata
f = fopen('nips17/nips17.docnames', 'w');
f2 = fopen('nips17/nips17.citations', 'w');
fc = fopen('nips17/nips17.labels', 'w');
fv = fopen('nips17/nips17.vols', 'w');
fy = fopen('nips17/nips17.years', 'w');
for i = 1:length(docs_names)
    % extract name year and volume number
    name = docs_names{i};
    year = eval(name(1:4));
    vol = year - 1987;
    cat = getcategory(name(6:end));
    % name
    fprintf(f, '%s\n', name);
    % more complete citation information
    fprintf(f2, '%04d %s: ', i - 1, name);
    aa = find(docs_authors(i, :));
    for j = 1:length(aa)
        fprintf(f2, '%s, ', authors_names{aa(j)});
    end
    fprintf(f2, 'NIPS %d, %d\n', vol, year);
    %fprintf(f2, 'end\n', vol, year);
    % remaining metadata
    fprintf(fy, '%d\n', year);
    fprintf(fv, '%d\n', vol);
    fprintf(fc, '%s\n', cat);
end
fclose(f);
fclose(f2);
fclose(fc);
fclose(fv);
fclose(fy);

% document authors
f = fopen('nips17/nips17.authors', 'w');
for i = 1:size(docs_authors, 1)
    % ndocs x nauth sparse, convert to zero offset
    aa = find(docs_authors(i, :)) - 1;
    for j = 1:length(aa)
        fprintf(f, '%d ', aa(j));
    end
    fprintf(f, '\n');
end
fclose(f);

% author names
f = fopen('nips17/nips17.authors.key', 'w');
for i = 1:length(authors_names)
    fprintf(f, '%s\n', authors_names{i});
end
fclose(f);

end %createnips

function cat = getcategory(string)
    c = string(1);
    d = string(3);
    % category pattern CC##, so check for alpha at 1 and num at 3
    if ( '0' <= d ) && ( d <= '9' ) && (('0' > c) || ('9' < c))
        cat = string(1:2);
    else
        cat = '--';
    end
end
