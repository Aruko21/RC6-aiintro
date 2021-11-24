% Представление базы фактов подмножества УСП

plate(30). % Высота плиты (мм)

% Таблица 1. Параметры опоры сферической (тип опоры, номенклатурный код, диаметр метр. резьбы)
support(sphere, '2.910.01', 6).
support(sphere, '2.910.02', 8).
support(sphere, '3.910.01', 12).
support(sphere, '3.910.02', 16).

% Таблица 2. Параметры опоры с насечкой (тип опоры, номенклатурный код, диаметр метр. резьбы (Dm))
support(notch, '2.911.01', 6).
support(notch, '2.911.02', 8).
support(notch, '3.911.01', 12).

% Таблица 3. Параметры штыря (номенклатурный код, диаметр метр. резьбы (Dm),  мин. диаметр (Dmin), макс. диаметр(Dmax))
pin('2.213.01', 6, 6, 8).
pin('2.213.04', 8, 8, 12).
pin('3.213.06', 12, 12, 26).

% Таблица 4. Параметры кулачка с резьбовым отверстием (тип кулачка, номенклатурный код, типоразмер, диаметр метр. резьбы (Dm), высота кулачка (H_кл))
tappet(carving, '2.913.05', '30_18', 6, 16).
tappet(carving, '2.913.06', '45_22', 8, 20).
tappet(carving, '2.913.09', '65_30', 12, 38).

% Таблица 5. Параметры кулачка с призматическими пазами (тип кулачка, номенклатурный код, типоразмер, высота кулачка (H_кл), диаметр горизонтальный мин/макс (Dhor_min/_max), диаметр вертикальный мин/макс (Dvert_min/_max))
tappet(prizm, '2.913.01', '30_18', 10, 8, 12, 3, 7).
tappet(prizm, '2.913.02', '45_22', 12, 8, 12, 3, 7).
tappet(prizm, '2.913.07', '65_30', 25, 12, 30, 8, 18).

% Таблица 6. Параметры зажима кулачкового (номенклатурный код, типоразмер зажима, типоразмер кулачка, высота зажима (H_зж))
clamp('2.451.01', '45_30', '30_18', 29).
clamp('2.451.02', '60_45', '45_22', 34).
clamp('3.451.01', '60_45', '45_22', 35).
clamp('3.451.02', '90_60', '65_30', 42).

% Таблица 7. Параметры прокладки (номенкенклатурный код, типоразмер, высота прокладки (H_пр)).
% Расположение прокладок по убыванию с целью оптимизации
gasket('3.107.28', '90_60', 5).
gasket('3.107.27', '90_60', 3).
gasket('3.107.25', '90_60', 2).

gasket('3.217.10', '60_45', 5).
gasket('3.217.09', '60_45', 3).
gasket('3.217.07', '60_45', 2).
gasket('3.217.01', '60_45', 1).

gasket('2.217.10', '45_30', 5).
gasket('2.217.09', '45_30', 3).
gasket('2.217.07', '45_30', 2).
gasket('2.217.01', '45_30', 1).

% --- Правила для проектирования пакета прокладок
gaskets_package([], _, Height, 0) :- Height = 0, !. % Требуемая высота = 0
gaskets_package([], _, Height, 0) :- Height < 0, !, fail. % Требуемая высота отрицательна

% Добавление прокладок в соотвтетсвии с типоразмером
gaskets_package([Code | GasketsListTail], TypeSize, PackageHeight, PackageAmount):-
	gasket(Code, TypeSize, GasketHeight),
	ReducedPackageHeight is PackageHeight - GasketHeight,
	gaskets_package(GasketsListTail, TypeSize, ReducedPackageHeight, ReducedPackageAmount),
	PackageAmount is ReducedPackageAmount + 1,
	!.

% --- Правила проектирования подвижной части приспособления

% Плоская черновая заготовка - опора с насечкой + кулачок с резьбовым отверстием (диаметр не требуется)
moving_part(flat_rough, _, CodeSupport, CodeTappet, TappetTypeSize, H_tappet) :-
	support(notch, CodeSupport, D), 
	tappet(carving, CodeTappet, TappetTypeSize, D, H_tappet).

% Плоская чистовая заготовка - сферическая опора + кулачок с резьбовым отверстием (диаметр не требуется)
moving_part(flat_clean, _, CodeSupport, CodeTappet, TappetTypeSize, H_tappet) :-
	support(sphere, CodeSupport, D), 
	tappet(carving, CodeTappet, TappetTypeSize, D, H_tappet).

% Цилиндрическая вертикальная заготовка - кулачок с призматическими пазами (опора не требуется)
moving_part(cyl_vert, D, _, CodeTappet, TappetTypeSize, H_tappet) :-
	nonvar(D), % Данная проверка необходима по причине невозможности проверки того, подходит ли заготовка, если переменная не конкретизированна
	tappet_prizm(CodeTappet, TappetTypeSize, H_tappet, _, _, D_min, D_max),
	D_max >= D, D_min =< D.

% Цилиндрическая горизонтальная заготовка - кулачок с призматическими пазами (опора не требуется)
moving_part(cyl_hor, D, _, CodeTappet, TappetTypeSize, H_tappet) :-
	nonvar(D), % Данная проверка необходима по причине невозможности проверки того, подходит ли заготовка, если переменная не конкретизированна
	tappet(prizm, CodeTappet, TappetTypeSize, H_tappet, D_min, D_max, _, _),
	D_max >= D, D_min =< D.

% Перфораторная заготовка - штырь + кулачок с резьбовым отверстием
moving_part(perf, D, CodeSupport, CodeTappet, TappetTypeSize, H_tappet) :-
	nonvar(D), % Данная проверка необходима по причине невозможности проверки того, подходит ли заготовка, если переменная не конкретизированна
	pin(CodeSupport, D, D_min, D_max),
	tappet_carving(CodeTappet, TappetTypeSize, D, H_tappet),
	D_max >= D, D_min =< D.


% --- Основной предикат, формирующий прижимное приспособление
% (Высота, )
device(SurfaceType, ClampHeight, Diameter) :-
    plate(PlateHeight),
	moving_part(SurfaceType, Diameter, CodeSupport, CodeTappet, TappetTypeSize, TappetHeight),
	clamp(CodeClamp, ClampTypeSize, TappetTypeSize, ClampHeightInternal),
	RequiredHeight is TappetHeight + ClampHeightInternal + PlateHeight,
	ResultHeight is ClampHeight - RequiredHeight,
	(
	    ResultHeight < 0 ->
	    write("Величины высоты зажима недостаточно.\nМинимальная требуемая высота для текущей конфигурации: "),
		write(RequiredHeight),
		nl, ! ; true
    ),
	gaskets_package(GasketCodeList, ClampTypeSize, RequiredHeight, GasketsAmount),
	print_device_info(SurfaceType, CodeSupport, CodeTappet, CodeClamp, RequiredHeight, GasketsAmount, GasketCodeList).

% Правила для вывода информации о спроектированном приспособлении
print_device_info(SurfaceType, CodeSupport, CodeTappet, CodeClamp, RequiredHeight, GasketsAmount, GasketCodeList) :-
	print_movable_part(SurfaceType, CodeSupport),
	print_tappet(SurfaceType, CodeTappet),
	print_clamp(CodeClamp),
	print_gaskets(RequiredHeight, GasketsAmount, GasketCodeList).

print_movable_part(flat_rough, CodeSupport) :-
	write("[Установочный элемент] - 'Опора с насечкой' | код: "), write(CodeSupport), nl.
print_movable_part(flat_clean, CodeSupport) :-
	write("[Установочный элемент] - 'Сферическая опора' | код: "), write(CodeSupport), nl.
print_movable_part(cyl_vert, _) :-
	write("[Установочный элемент] - не требуется"), nl.
print_movable_part(cyl_hor, _) :-
	write("[Установочный элемент] - не требуется"), nl.
print_movable_part(perf, CodeSupport) :-
	write("[Установочный элемент] - 'Штырь' | код: "), write(CodeSupport), nl.

print_tappet(flat_rough, CodeTappet) :-
	write("[Направляющий элемент] - 'Кулачок с резьбовым отверстием' | код: "), write(CodeTappet), nl.
print_tappet(flat_clean, CodeTappet) :-
	write("[Направляющий элемент] - 'Кулачок с резьбовым отверстием' | код: "), write(CodeTappet), nl.
print_tappet(cyl_vert, CodeTappet) :-
	write("[Направляющий элемент] - 'Кулачок с призматическими пазами' | код: "), write(CodeTappet), nl.
print_tappet(cyl_hor, CodeTappet) :-
	write("[Направляющий элемент] - 'Кулачок с призматическими пазами' | код: "), write(CodeTappet), nl.
print_tappet(perf, CodeTappet) :-
	write("[Направляющий элемент] - 'Кулачок с резьбовым отверстием' | код: "), write(CodeTappet), nl.

print_clamp(CodeClamp) :-
	write("[Зажимной элемент] - 'Зажим' | код: "), write(CodeClamp), nl.

print_gaskets(RequiredHeight, GasketsAmount, GasketCodeList) :-
    write("[Пакет прокладок]:"), nl,
	write("Высота: "), write(RequiredHeight), nl,
	write("Общее количество: "), write(GasketsAmount), nl,
	write("Элементы:"), nl,
	print_gaskets_list(_, GasketCodeList).

print_gaskets_list(OutList, GasketsList) :-
	uniqlist(GasketsList, UniqueGaskets),
	gasket_groups(AmountList, GasketsList, UniqueGaskets),
	print_gasket_groups(OutList, AmountList, UniqueGaskets),
	!.		

uniqlist([], []) :- !.
uniqlist([H | T], Out) :- member(H, T), uniqlist(T, Out).
uniqlist([H | T], [H | T1]) :- not(member(H, T)), uniqlist(T, T1).

% Формирование групп подкладок с подсчетом элементов в каждой
gasket_groups([], _, []).
gasket_groups([GroupCount | GroupCountTail], GasketsList, [HeadUniqGasket | UniqueGasketsTail]) :-
	add_gaskets_to_group(GroupCount, GasketsList, HeadUniqGasket),
	gasket_groups(GroupCountTail, GasketsList, UniqueGasketsTail).

add_gaskets_to_group(0, [], _).
add_gaskets_to_group(GroupCount, [Gasket | GasketsTail], GasketGroup) :-
	GasketGroup \= Gasket, %( != )
	add_gaskets_to_group(GroupCount, GasketsTail, GasketGroup).
add_gaskets_to_group(GroupCount, [Gasket | GasketsTail], GasketGroup) :-
	GasketGroup = Gasket,
	add_gaskets_to_group(GroupCountPrev, GasketsTail, GasketGroup),
	GroupCount is GroupCountPrev + 1.

print_gasket_groups([], [], []).
print_gasket_groups([[UniqueGasket : Amount] | GasketsTail], [Amount | AmountsTail], [UniqueGasket | UniqueGasketsTail]) :-
	write("Код: "), write(UniqueGasket), write("; количество: "), write(Amount), nl,
	print_gasket_groups(GasketsTail, AmountsTail, UniqueGasketsTail).
