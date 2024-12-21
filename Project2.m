clear clc;
%Part A
disp("Hello welcome to Santiago's: Project 2 code ");
disp(' ');
file=input("PART A: \n Please enter a valid txt file in the command window : \n ","s");
disp(' ');
input_file =fopen(file,'r'); % read the input file
%input_file=strtrim(fscanf(file));
if(input_file == -1)
    error("File can't be read please try checking your spelling or making sure the file is in the current folder");
end

file_strand=strtrim(fscanf(input_file,'%c'));
fclose(input_file);
%Part B
choice=input("PART B: Thank you for entering the file \nIs the file you provided a coding or template strand? Please enter either (coding or template: )","s");
disp(' ');
randomval=0;
choice1="I see you are outputting a "+choice+" strand";
disp(choice1);
choice_of_direction=input("Please let me know what direction the strand is going \neither type 5-3 or 3-5\n","s");
choice_of_direction1=choice_of_direction;
if strcmpi(choice, "coding")
    randomval=1;
    choice_of_direction1='3-5';
    if strcmpi(choice_of_direction, '3-5')
     file_strand=reverse(file_strand);
     choice_of_direction1='3-5';
    end
end
if strcmpi(choice,"template")
     choice_of_direction1='5-3';
    if strcmpi(choice_of_direction, '5-3')
     file_strand=reverse(file_strand);
     choice_of_direction1='5-3';
    end
end

    

%Part C
disp("PART C: Thank you for inputing the information, I will now display the complentary strand for your given file:")
disp(' ');
output_strand_name="Complemtary strand in "+choice_of_direction1+" orientation:";
%disp("Complentary strand:");
disp(output_strand_name);
    complement_line=file_strand;
    complement_line = regexprep(complement_line, '[^ATGC]', '');
    complement_line=strrep(complement_line,'A','1');
    complement_line=strrep(complement_line,'T','2');
    complement_line=strrep(complement_line,'G','3');
    complement_line=strrep(complement_line,'C','4');
    complement_line=strrep(complement_line,'1','T');
    complement_line=strrep(complement_line,'2','A');
    complement_line=strrep(complement_line,'3','C');
    complement_line=strrep(complement_line,'4','G');
    fprintf("%s\n", complement_line);
    
    
%part D
    %disp("Debugging: Checking complement_line");
    %disp(complement_line);
    %disp("Length of complement_line: " + length(complement_line));
    disp("PART D: The number of each nucleic acid in each complimentary strand");
    disp(' ');
    if isempty(complement_line)
        error("Please use a txt file with valid nucleotides")
    end

    number_of_A=sum(complement_line== 'A');
    number_of_T=sum(complement_line== 'T');
    number_of_C=sum(complement_line== 'C');
    number_of_G=sum(complement_line== 'G');
    fprintf("A: %d\nT: %d\nC: %d\nG: %d\n", number_of_A, number_of_T, number_of_C, number_of_G);
    counts = [number_of_A, number_of_T, number_of_C, number_of_G];
    nucleotide={'A','T','C','G'};
    
    figure;    
    bar(counts, 'FaceColor', 'blue');
    set(gca, 'XTickLabel', nucleotide, 'XTick', 1:numel(nucleotide));
    xlabel('Nucleotide');
    ylabel('Count');
    title('Distribution of Nucleotides for Complementary Strand');
 %part E
    number_of_T_comp=sum(complement_line== 'A');
    number_of_A_comp=sum(complement_line== 'T');
    number_of_G_comp=sum(complement_line== 'C');
    number_of_C_comp=sum(complement_line== 'G');

    num_total_A=number_of_A+number_of_A_comp;
    num_total_T=number_of_T+number_of_T_comp;
    num_total_C=number_of_C+number_of_C_comp;
    num_total_G=number_of_G+number_of_G_comp;
    disp("PARTE: Total number of nucleotides in dsDNA:")
    fprintf("A: %d\nT: %d\nC: %d\nG: %d\n", num_total_A, num_total_T, num_total_C, num_total_G);
    counts1 = [num_total_A, num_total_T, num_total_C, num_total_G];
    nucleotide1={'A','T','C','G'};
    figure;    
    bar(counts, 'FaceColor', 'blue');
    set(gca, 'XTickLabel', nucleotide1, 'XTick', 1:numel(nucleotide1));
    xlabel('Nucleotide');
    ylabel('Count');
    title('Distribution of Nucleotides for dsDNA Strand');
    %title('Distribution of Nucleotides for dsDNA');
%part F
    pairs_of_AT=number_of_T*2;
    pairs_of_GC=number_of_C*3;
    hydrogen_bonds=pairs_of_GC+pairs_of_AT;
    disp(['PART F:Hydrogren bonds: ', num2str(hydrogen_bonds)]);
    %disp(['Pairs Of AT: ', num2str(pairs_of_AT)]);
    %disp(['Pairs Of GC: ', num2str(pairs_of_GC)]);
%part G
%equation for MW
%M.W. of dsDNA = (# nucleotides x 607.4) + 157.9
%reference1 https://www.thermofisher.com/us/en/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html

MW_dsDNA=((num_total_A+num_total_T+num_total_G+num_total_C)* 607.4 + 157.9)/2;
disp(['PART G: MW of dsDNA: ', num2str(MW_dsDNA)]);
%part H
    if randomval==1
        mRNA=complement_line;
        mRNA = regexprep(mRNA, '[^ATGC]', '');
        mRNA=strrep(mRNA,'A','1');
        mRNA=strrep(mRNA,'T','2');
        mRNA=strrep(mRNA,'C','3');
        mRNA=strrep(mRNA,'G','4');

        mRNA=strrep(mRNA,'1','U');
        mRNA=strrep(mRNA,'2','A');
        mRNA=strrep(mRNA,'3','G');
        mRNA=strrep(mRNA,'4','C');

        %mRNA=strrep(mRNA,'2','G');
    end
    if randomval==0
        %file_strand=reverse(file_strand);
        mRNA=file_strand;
        mRNA = regexprep(mRNA, '[^ATGC]', '');
        mRNA=strrep(mRNA,'A','1');
        mRNA=strrep(mRNA,'T','2');
        mRNA=strrep(mRNA,'C','3');
        mRNA=strrep(mRNA,'G','4');

        mRNA=strrep(mRNA,'1','U');
        mRNA=strrep(mRNA,'2','A');
        mRNA=strrep(mRNA,'3','G');
        mRNA=strrep(mRNA,'4','C');
    end
    disp(' ');
    disp('PART H: mRNA of coding strand: ')
    fprintf("%s\n", mRNA);
    %part I
    map_of_codons=containers.Map({ ...
    'UUU', 'UUC', 'UUA', 'UUG', ...
    'CUU', 'CUC', 'CUA', 'CUG', ...
    'AUU', 'AUC', 'AUA', 'AUG', ...
    'GUU', 'GUC', 'GUA', 'GUG', ...
    'UCU', 'UCC', 'UCA', 'UCG', ...
    'CCU', 'CCC', 'CCA', 'CCG', ...
    'ACU', 'ACC', 'ACA', 'ACG', ...
    'GCU', 'GCC', 'GCA', 'GCG', ...
    'UAU', 'UAC', 'UAA', 'UAG', ...
    'CAU', 'CAC', 'CAA', 'CAG', ...
    'AAU', 'AAC', 'AAA', 'AAG', ...
    'GAU', 'GAC', 'GAA', 'GAG', ...
    'UGU', 'UGC', 'UGA', 'UGG', ...
    'CGU', 'CGC', 'CGA', 'CGG', ...
    'AGU', 'AGC', 'AGA', 'AGG', ...
    'GGU', 'GGC', 'GGA', 'GGG'}, ...
    {'PHE','PHE','LEU','LEU' ...
    ,'LEU','LEU','LEU','LEU' ...
    ,'LLE','LLE','LLE','MET' ...
    ,'VAL','VAL','VAL','VAL' ...
    ,'SER','SER','SER','SER' ...
    ,'PRO','PRO','PRO','PRO' ...
    ,'THR','THR','THR','THR' ...
    ,'ALA','ALA','ALA','ALA' ...
    ,'TYR','TYR','STOP','STOP' ...
    ,'HIS','HIS','GLN','GLN' ...
    ,'ASN','ASN','LYS','LYS' ...
    ,'ASP','ASP','GLU','GLU' ...
    ,'CYS','CYS','STOP','TRP' ...
    ,'ARG','ARG','ARG','ARG' ...
    ,'SER','SER','ARG','ARG' ...
    ,'GLY','GLY','GLY','GLY'});
amino_acid_counter = containers.Map('KeyType', 'char', 'ValueType', 'double');
amino_acid_chain = " ";
for i = 1:3:length(mRNA)
    codons = mRNA(i:min(i+2, length(mRNA)));
   if ~isKey(map_of_codons, codons)
        fprintf("Invalid codon encountered:", codons); %#ok<*CTPCT>
        continue; % Skip invalid codons
    end
    amino_acid = map_of_codons(codons);
    disp(amino_acid);
    amino_acid_chain = amino_acid_chain +" "+ amino_acid;
    %if amino_acid == "STOP"
        %amino_acid_chain = amino_acid_chain +" "+ amino_acid;
        
    %end
    if amino_acid_counter.isKey(amino_acid)
        amino_acid_counter(amino_acid) = amino_acid_counter(amino_acid) + 1;
    else
        amino_acid_counter(amino_acid) = 1; % Initialize counter if first occurrence
    end
end
disp('PART I:')
disp(amino_acid_chain);
disp('PART J')
disp('Amino acid counts:');
keys = amino_acid_counter.keys;
values = amino_acid_counter.values;
for j = 1:length(keys)
    fprintf('%s: %d\n', keys{j}, values{j});
end
 figure;
    bar(cell2mat(values), 'FaceColor', 'red');
    set(gca, 'XTickLabel', keys, 'XTick', 1:numel(keys));
    xlabel('Amino Acid');
    ylabel('Count');
    title('Distribution of Amino Acids');