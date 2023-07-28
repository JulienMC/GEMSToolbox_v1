function [tf] = input_eval_check(input)
    valid_chars=["0","1","2","3","4","5","6","7","8","9",".",";","[","]","_","e","E","-"];  % allowed characters
    allCharsGood = @(c, goodList)all(ismember(string(cellstr(c')'), goodList)); 
    input = strrep(input,' ','');
    tf = allCharsGood(input,valid_chars);
end