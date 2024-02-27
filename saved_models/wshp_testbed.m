% Example input data
Xnew = [0.5, 0.2, 15, 25, 50;  % Example row 1
        0.3, 0.1, 16, 26, 55]; % Example row 2

% Load the model using Python command
model = py.joblib.load('wshp_model.joblib');

% Initialize a cell array to store predictions
y_pred = cell(1, size(Xnew, 1));

% Convert each row of Xnew to a Python list and predict
for i = 1:size(Xnew, 1)
    % Extract a single row and convert it to a MATLAB list
    row = Xnew(i, :);
    
    % Convert the MATLAB row to a Python list
    py_row = py.list(row);
    
    % Convert the MATLAB list to a nested Python list as expected by the predict method
    py_nested_list = py.list({py_row});
    
    % Predict using the model and store the result
    y_pred{i} = model.predict(py_nested_list);
end

% Display the prediction results
disp(y_pred);
