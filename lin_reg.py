# ---- Åsmunds kokte funksjon for lineær regresjon ---- 
# Input er en x- og y-vektor med samme dimensjon
# Output er en linearisert y-vektor som kan plottes som funksjon av den opprinnelige x-vektor
# Printer også funksjonsuttrykket

import numpy as np

def lin_reg(x,y):
    x=x.reshape((-1, 1))
    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    model.fit(x, y)
    model = LinearRegression().fit(x, y)
    r_sq = model.score(x, y)
    y_pred = model.predict(x)
    y_pred = model.intercept_ + np.sum(model.coef_ * x, axis=1)
    #x_new = np.arange(10).reshape((-1, 2))
    linear_y = model.predict(x)


    stigtall = (linear_y[-1]-linear_y[0])/(x[-1]-x[0])
    konstledd = linear_y[-1]-stigtall*x[-1]
    print('Funksjonsuttrykket blir y=',stigtall,'x +',konstledd)
    return linear_y