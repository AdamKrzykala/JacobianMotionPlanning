## Jacobian motion planning algorithm

Other types of jacobian inverse, for instance:
- pseudoinverse (Moore - Penrose)
- Extended lagrangian jacobian

Program ECSA 2:

1. W celu przeprowadzenia zadania planowania ruchu należy edytować tylko i wyłącznie plik 
taskDef.py. 
2. Edytowalne zmienne w programie:
	- Th: horyzont czasowy ruchu 
	- ConfDim: wymiar wektora zmiennych stanu
	- ControlDim: ilość sterowań 
	- OutputDim: wymiar wekotra wyjścia systemu
	- controlVectorSize: wymiar wekotra parametrów funkcji sterujących. 
		(jeden parametr określa zawsze wyraz wolny, dwa kolejne pierwsze harmoniczne, 
		dwa ostatnie drugie harmoniczne). W przykładzie podano 10, co oznacza, że są dwa 
		sterowania, każde ma 5 parametrów - wyraz wolny, pierwsze i drugie harmoniczne). 
		Przy dwóch sterowaniach możliwe wartości:
		2 - tylko wyrazy wolne, 
		6 - wyraz wolny + pierwsze harmoniczne
		10 - wyraz wolny + pierwsze harmoniczne + drugie harmoniczne
	- InitConfiguration: wartości początkowe zmiennych stanów na początku ruchu
	- InitControl: wartości początkowe wektora parametrów funkcji sterujących
	- Część USER VARIABLES pozwala na zdefiniowanie zmiennych użytkownika
	- Generators Matrix - macierz G(q)
	- initialDrift - dryf układu f(q)
	- outputSys - przekształcenie zmiennych stanu w wyjście 
	- InitGamma - początkowa wartość parametru gamma 
	- desPos - zadane wyjście systemu, które musi być osiągnięte po czasie Th
	- maxError - maksymalny dopuszczalny błąd

3. Uruchamianie programu:
	python3.6 main.py

4. Założenia programu: 
	- Odwrotnośc lagranżowska (możliwość implementacji różnych odwrotności)
	- Sterowanie na początku i końcu osiąga wartość 0
	- funkcje sterujące występują parami (sin + cos) - opcja używania pojedynczo
	
5. Efekt działania programu:
	- Step: numer kroku
	- Control value at the beginning - powinno być bliskie zera - założone gładkie sterowania 
	- Control value at the end - powinno być bliskie zera - założone gładkie sterowania 
	- Configuration: wektor konfiguracji w danym kroku
	- Error: błąd pomiędzy wyjściem aktualnym, a założonym

6. Dalszy rozwój:
	- Zwiększenie elastyczności programu, 
	- Interfejs graficzny

Author: Adam Krzykala
