/* Задание 9. Вариант 50
	Написать программу символьного дифференцирования арифметических выражений, представленных в канонической (в виде структур)
	для языка Пролог форме. При этом единственная переменная (по которой производится дифференцирование) обозначается буквой x,
	а константы - всеми иными строчными (маленькими) буквами.
*/

% Система символьного диффренцирования представляет собой набор правил для каждой атомарной операции в отдельности

% Операция сложения
simp_sum(A, 0, A).
simp_sum(0, A, A).

simp_sum(A, B, Sum) :-
	number(A),
	number(B),
	Sum is A + B.

simp_sum(A, B, A + B).

% Операция вычитания
simp_sub(A, 0, A).
simp_sub(0, A, -A).

simp_sub(A, B, Sub) :-
	number(A),
	number(B),
	Sub is A - B.

simp_sub(A, B, 0) :- A == B.

simp_sub(A, B, A - B).

% Операция умножения
simp_mul(0, _, 0).
simp_mul(_, 0, 0).

simp_mul(A, 1, A).
simp_mul(1, A, A).

simp_mul(A, B, Mult) :-
	number(A),
	number(B),
	Mult is A * B.

simp_mul(A, B, A * B).

% Операции деления
simp_div(_, 0, _) :- throw("Divide by zero error").

simp_div(_, 0, 0).

simp_div(A, 1, A).

simp_div(A, A, 1).

simp_div(A, B, Div) :-
	number(A),
	number(B),
	Div is A / B.

simp_div(A, B, A / B).


% Операция дифференцирования

% d(const)/dx = 0
diff(Eq, _, 0) :- number(Eq).

% d(y)/dx = 0
diff(Eq, Var, 0) :- 
	atom(Eq),
	Eq \== Var.

% d(x)/dx = 1
diff(Eq, Var, 1) :- 
	atom(Eq),
	Eq == Var.

% d(u + v)/dx = d(u)/dx + d(v)/dx
diff(+(Eq1, Eq2), Var, Res) :-
	diff(Eq1, Var, Res1),
	diff(Eq2, Var, Res2),
	simp_sum(Res1, Res2, Res).

% d(u - v)/dx = d(u)/dx - d(v)/dx
diff(-(Eq1, Eq2), Var, Res) :-
	diff(Eq1, Var, Res1),
	diff(Eq2, Var, Res2),
	simp_sub(Res1, Res2, Res).

% d(u * v)/dx = v * (d(u)/dx) + u * (d(v)/dx)
diff(*(Eq1, Eq2), Var, Res) :-
	diff(Eq1, Var, Res1),
	diff(Eq2, Var, Res2),
	simp_mul(Eq2, Res1, L),
	simp_mul(Eq1, Res2, R),
	simp_sum(L, R, Res).

% d(u / v)/dx = (v * (d(u)/dx) - u * (d(v)/dx)) / (v * v)
diff(/(Eq1, Eq2), Var, Res) :-
	diff(Eq1, Var, Res1),
	diff(Eq2, Var, Res2),
	simp_mul(Eq2, Res1, L),
	simp_mul(Eq1, Res2, R),
	simp_mul(Eq2, Eq2, Denom),
	simp_sub(L, R, Nom),
	simp_div(Nom, Denom, Res).

% d(sin(u))/dx = cos(u) * (d(u)/dx)
diff(sin(Eq), Var, Res) :-
    diff(Eq, Var, Res1),
	simp_mul(cos(Eq), Res1, Res).

% d(cos(u))/dx = -sin(u) * (d(u)/dx)
diff(cos(Eq), Var, Res) :-
    diff(Eq, Var, Res1),
	simp_mul(sin(Eq), Res1, Mult),
	simp_sub(0, Mult, Res).

% d(tg(u))/dx = (1 / cos(u)^2) * (d(u)/dx) 
diff(tg(Eq), Var, Res) :-
	diff(Eq, Var, Res1),
	simp_div(Res1, cos(Eq)^2, Res).

% d(exp(u))/dx = exp(u) * (d(u)/dx)
diff(exp(Eq), Var, Res) :-
    diff(Eq, Var, Res1),
	simp_mul(exp(Eq), Res1, Res).

% d(ln(u))/dx = (1 / u) * (d(u)/dx)
diff(ln(Eq), Var, Res) :-
    diff(Eq, Var, Res1),
	simp_div(Res1, Eq, Res).

% d(f(u)^a)/dx = a * f(u)^(a - 1) * (d(f(u))/dx)
diff(^(Eq, Num), Var, Res) :-
    number(Num),
    diff(Eq, Var, Res1),
    simp_sub(Num, 1, Pow),
    simp_mul(Num, (Eq)^Pow, R),
    simp_mul(R, Res1, Res).

% d(a^(f(u)))/dx = a^(f(u)) * ln(a) * (d(f(u))/dx)
diff(^(Num, Eq), Var, Res) :-
    number(Num),
    diff(Eq, Var, Res1),
    simp_mul(Num^(Eq), ln(Num), R),
    simp_mul(R, Res1, Res).

% d(f(u)^g(u))/dx = f(x)^(g(x) - 1) * (g(x) * d(f(x))/dx + d(g(x))/dx * f(x) * ln(f(x)) 
diff(^(Eq1, Eq2), Var, Res) :-
    diff(Eq1, Var, Res1),
    diff(Eq2, Var, Res2),
	simp_sub(Eq2, 1, Pow),
	simp_mul(Eq2, Res1, L),
	simp_mul(Eq1, ln(Eq1), Tmp),
	simp_mul(Tmp, Res2, R),
	simp_sum(L, R, Sum),
	simp_mul(Eq1^(Pow), Sum, Res).
