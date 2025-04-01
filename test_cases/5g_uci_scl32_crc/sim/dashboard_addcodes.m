result_files = dir('results_polar*L=32*.txt')
code_files = dir('code_polar*L=32*.txt')

for ii=1:length(code_files)
    fprintf('python3 data_import.py -r ../codes/polar/5g_uci_scl32_crc/sim/%s -d SCL -c ../codes/polar/5g_uci_scl32_crc/sim/%s --cat 18  --mod ASK --modset 2 --addedby thomas_wiegart\n', result_files(ii).name, code_files(ii).name);
end