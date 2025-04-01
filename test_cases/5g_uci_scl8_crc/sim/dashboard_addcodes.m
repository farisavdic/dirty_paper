result_files = dir('results_polar*L=8_*.txt')
code_files = dir('code_polar*L=8*.txt')

for ii=1:length(result_files)
    fprintf('python3 data_import.py -r ../codes/polar/5g_uci_scl8_crc/sim/%s -d SCL -c ../codes/polar/5g_uci_scl8_crc/sim/%s --cat 17  --mod ASK --modset 2 --addedby thomas_wiegart\n', result_files(ii).name, code_files(ii).name);
end