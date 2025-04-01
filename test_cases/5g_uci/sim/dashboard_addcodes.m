result_files = dir('results_polar*L=1_CRC-0.txt')
code_files = dir('code_polar*L=1_CRC-0.txt')

for ii=1:length(result_files)
    fprintf('python3 data_import.py -r ../codes/polar/5g_uci/sim/%s -d SC -c ../codes/polar/5g_uci/sim/%s --cat 16  --mod ASK --modset 2 --addedby thomas_wiegart\n', result_files(ii).name, code_files(ii).name);
end