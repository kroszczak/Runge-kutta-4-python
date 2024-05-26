def infix_to_rpn(expression):
    # Słownik priorytetów operatorów
    precedence = {'+': 1, '-': 1, '*': 2, '/': 2, '^': 3}

    output = []  # Wynikowa notacja RPN
    stack = []   # Stos operatorów

    # Funkcja pomocnicza do obsługi operatorów
    def handle_operator(operator):
        while stack and precedence.get(operator, 0) <= precedence.get(stack[-1], 0):
            output.append(stack.pop())
        stack.append(operator)

    # Przechodzenie po każdym tokenie w wyrażeniu infixowym
    for token in expression:
        if token.isdigit():
            output.append(token)
        elif token in precedence:
            if token == '^':
                # Obsługa potęgowania jako operatora z najwyższym priorytetem
                stack.append(token)
            else:
                handle_operator(token)
        elif token == '(':
            stack.append(token)
        elif token == ')':
            while stack and stack[-1] != '(':
                output.append(stack.pop())
            stack.pop()  # Usuń '(' ze stosu
            if stack and stack[-1] == '^':
                output.append(stack.pop())  # Dodaj operator potęgowania na wyjście
    while stack:
        output.append(stack.pop())
    input(output)
    return ' '.join(output).replace("^", "**")