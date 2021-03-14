%jjs280
%3/14/2021
%Reads variables from an input file

function Inputs = read_inputs(Executor)
    %Specify input file
    disp(Executor);
    inputFile = input('enter input file -> ','s');
    try
        fin = fopen(inputFile,'r');
    catch
        disp('Could not open input file');
        disp(' ');
        Inputs = 0;
        return
    end
    disp(' ');

    % Read inputs from file and stores them in a struct
    s = fgetl(fin);
    
    
end