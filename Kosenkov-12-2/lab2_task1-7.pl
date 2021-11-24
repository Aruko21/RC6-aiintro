% Задание 1 - Конкатенация списков: concat(L1,L2,L3).
% Предикат удовлетворен, если список L3 является конкатенацией (слиянием) списков L1 и L2

concat([], ListConcat, ListConcat).
concat([Head1 | List1Tail], List2, [Head1 | ConcatTail]) :- concat(List1Tail, List2, ConcatTail).

% Задание 2 - Инверсия списков: invert(L1,L2).
% Предикат удовлетворен, если список L2(L1) является инверсией списка L1(L2).

invert(ListInit, ListReverse) :- check_reverse(ListInit, [], ListReverse).
check_reverse([], ListReverse, ListReverse).
check_reverse([HeadInit | ListInitTail], ListReverse, ListCheck) :-
    check_reverse(ListInitTail, [HeadInit | ListReverse], ListCheck).

% Задание 3 - Элементы списка без повторений: uniq(L1,L2).
% Предикат удовлетворен, если список L2 содержит все элементы L1 без повторений.

uniq(L1, L2) :- invert(L1, L1Invert), invert(L2, L2Invert), uniqlist(L1Invert, L2Invert).

uniqlist([], []) :- !.
uniqlist([H | T], Out) :- nonvar(H), member(H, T), uniqlist(T, Out).
uniqlist([HList | T], [HUniq | T1]) :- nonvar(HList), HList = HUniq, not(member(HList, T)), uniqlist(T, T1).
uniqlist([HList | T], [HUniq | T1]) :- var(HList), HList = HUniq, uniqlist(T, T1). % Данное правило необходимо для корректного формирования List1 в случае, если список неконкретизирован полностью

member(_, []) :- !, fail. % member - Проверка принандлежности элемента списку
member(X, [ListHead | _]) :- nonvar(ListHead), X = ListHead, !. % nonvar ипользован для защиты от случаев, когда список состоит из неконкретизированных переменных
member(X, [_ | T]) :- member(X, T). % Переход к следующему элементу

% Задание 4 - Конкатенация первого списка и элементов второго, не присутствующих в первом: ucat(L1,L2,L3).
% Предикат удовлетворен, если список L3 содержит все элементы списка L1 и те элементы из списка L2, которых нет в L1.

ucat(List1, List2, ListUniq) :- ucat_reverse(ListUniq, List2, List1), filter(List2, List1, List2Filter), concat(List1, List2Filter, ListUniq).

ucat_reverse(ListUniq, List2, List1) :- % Формирование первого списка, удовлетворяющего правилам предиката
    % Алгоритм - от результирующей конкатенации с конца "отрезаются" элементы 2-го списка, порядок которых соответсвует порядку элементов второго списка
    invert(ListUniq, ListUniqReverse), invert(List2, List2Reverse), invert(List1, List1Reverse),
    ucat_reverse_internal(ListUniqReverse, List2Reverse, List1Reverse).
ucat_reverse_internal(List1, [], List1).
ucat_reverse_internal(ListUniq, [_ | List2Tail], List1) :- % Если на конце результирующего списка нет элемента, идущего по порядку во втором списке - значит он присутствует в первом
    ucat_reverse_internal(ListUniq, List2Tail, List1).
ucat_reverse_internal([Head | ListUniqTail], [Head | FilterListTail], List1) :- % Иначе - он уникален для второго
    ucat_reverse_internal(ListUniqTail, FilterListTail, List1).

filter([], _, []) :- !. % Удаление из второго списка элементов, встречающихся в первом
filter([Head | ListTail], FilterList, [Head | ResultTail]) :- not(member(Head, FilterList)), filter(ListTail, FilterList, ResultTail).
filter([Head | ListTail], FilterList, ResultList) :- member(Head, FilterList), filter(ListTail, FilterList, ResultList).

% Задание 5 - Отображение двух списков операцией: mapop(Op,L1,L2,L3).
% Предикат удовлетворен, если i-ый элемент списка L3 представляет собой результат применения инфиксной операции Op к i-ым элементам списков L1 и L2.

mapop(_, [], [], []) :- !.
mapop(Op, [Head1 | List1Tail], [Head2 | List2Tail], [HeadResult | ListResultTail]) :-
    Functor =..[Op, Head1, Head2], HeadResult is Functor, mapop(Op, List1Tail, List2Tail, ListResultTail).

% Задание 6 - Преобразование к одноуровневому списку: unbr(L1,L2).
% Предикат удовлетворен, если L1 есть многоуровневый список произвольных элементов, а L2 - одноуровневый список тех же элементов.

unbr([], []) :- !.
unbr([Head1 | List1Tail], List2) :- !, unbr(Head1, FlatHead1), unbr(List1Tail, FlatTail1), concat(FlatHead1, FlatTail1, List2).
unbr(List, [List]).

% Задание 7 - Сумма элементов подсписков списка: msum(L1,L2).
% Предикат удовлетворен, если L1 есть список списков чисел (произвольной длины), а L2 - список сумм чисел во вложенных списках.

msum(List1, List2) :- msum(List1, [], List2).
msum([], Result, List2) :- invert(Result, List2), !.
msum([Head1 | List1Tail], Result, List2) :- sum_list(Head1, Sum), msum(List1Tail, [Sum | Result], List2).

sum_list([], 0).
sum_list([Head | ListTail], Sum) :- sum_list(ListTail, SumTmp), Sum is SumTmp + Head.
