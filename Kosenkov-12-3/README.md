# Лабораторная работа №3 по курсу "Введение в искусственный интеллект"

## Вариант 1-7
## 1 - персептрон
## 7-я обучающая выборка

## Файлы программы
- Makefile - Makefile для компиляции и вспомогательных целей
- Generator.hpp - заголовочный файл генератора, используемого для создания наборов точек
- Perceptron.hpp - заголовочный файл персептрона
- Generator.cpp - файл с исходным кодом генератора, используемого для создания наборов точек
- Perceptron.cpp - файл с исходным кодом персептрона
- main.cpp - основной исполняемый файл

## Команды для компиляции
`make` или `make lab3`

## Команда для запуска программы
`./lab3`

## После того, как программа отработала, можно построить анимацию процесса обучения:
- Анимация в отдельном окне, для перехода к последующему кадру нужно нажать кнопку мыши:
`make plot_window`

- Зацикленная анимация, сохраняемая в файл result.gif
`make plot_gif`

- На графике отображена линия, разбивающая пространство входных сигналов на два класса

- Квадратные точки обозначают данные обучающей выборки, точки в виде зведочки обозначают тестовые данные

- Красным цветом обозначаются точки, принадлежащие к одному классу, синим - к другому

- В левом нижнем углу отображается номер итерации весов (с учетом всех переобучений)