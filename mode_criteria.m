criteria= menu('Mode validation criteria','Modal Quality',' Uncertainty','DP','Exit');
switch criteria
    case 1
        MQ
    case 2
        Uncertainty_criteria
    case 3
        DirichletProcess

end


if criteria <4
    mode_criteria
end
