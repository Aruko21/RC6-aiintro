:- set_prolog_flag(answer_write_options, [max_depth(0)]).

/* 

Пример: ?-path([1, 1], [5, 7], Path).
Path = [[1,1],[2,1],[3,1],[3,2],[4,2],[4,3],[4,4],[4,5],[4,6],[4,7],[5,7]].

Карта лабиринта.
Условные обозначения:
. - комната
# - преграда

 X\Y 1 2 3 4 5 6 7
  1  . . . . # . .
  2  . # . . # . .
  3  . # # # . . .
  4  . . . . . # .
  5  . . . . . # .
  6  . # # # . # .
  7  . . . # . # .

*/

% Задание размера лабиринта
maze_map_size(7, 7).

% Задание расположения преград в лабиринте: wall(X, Y)
wall(1, 5).
wall(2, 2).
wall(2, 5).
wall(3, 2).
wall(3, 3).
wall(3, 4).
wall(4, 6).
wall(5, 6).
wall(6, 2).
wall(6, 3).
wall(6, 4).
wall(6, 6).
wall(7, 4).
wall(7, 6).

% Поиск пути
% From = [X1, Y1], To = [X2, Y2], Path = [From, ..., To]
path(From, To, [From | Path]) :-
    find_path(From, To, [], Path).    

% Предикат поиска пути, пока Next \= From
find_path(To, To, Path, PrintPath) :- reverse(Path, PrintPath), !.

find_path(From, To, Path, PrintPath) :-
    is_valid(From, Next),
    not(member(Next, Path)),
    find_path(Next, To, [Next | Path], PrintPath).    

% Проверка зацикливания
member(_, []) :- !, fail.
member(X, [X | _]) :- !.
member(X, [_ | T]) :- member(X, T).

% Проверка двух координат
is_valid(From, To) :-
    adjacent(From, To),
    in_bound(From),
    in_bound(To),
    not_wall(From),
    not_wall(To).

% Правила для проверки комнат на соседство
adjacent([X1, Y1], [X2, Y2]) :-
    adjacent(X1, Y1, X2, Y2).

adjacent(Y, X1, Y, X2) :-
    X2 is X1 + 1 ;
    X2 is X1 - 1.

adjacent(Y1, X, Y2, X) :-
    Y2 is Y1 + 1 ; 
    Y2 is Y1 - 1.
 
% Проверка на нахождение координат в границах лабиринта
in_bound([H, W]) :-
    maze_map_size(MaxHeight, MaxWidth),
    W =< MaxWidth,
    W > 0,    
    H =< MaxHeight,
    H > 0.

% Проверка на наличие барьера в координате
not_wall([H, W]) :-
    not(wall(H, W)).

is_valid_room(Room):-
  in_bound(Room),
  not_wall(Room).
