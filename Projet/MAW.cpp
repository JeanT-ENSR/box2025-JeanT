#include <iostream>
#include <seqan/index.h>

using namespace seqan2;

int main ()
{
    String<char> myString = "ACTAACTG";

    typedef Index< String<char> > TMyIndex;
    TMyIndex myIndex(myString);

    typedef Iterator< TMyIndex, MaxRepeats >::Type TMaxRepeatIterator;
    TMaxRepeatIterator myRepeatIterator(myIndex, 3);

    while (!atEnd(myRepeatIterator)) 
    {
        Iterator<TMaxRepeatIterator>::Type myRepeatPair(myRepeatIterator);
        while (!atEnd(myRepeatPair)) {
            ::std::cout << *myRepeatPair << ", ";
            ++myRepeatPair;
        }

        ::std::cout << repLength(myRepeatIterator) << "   ";

        ::std::cout << "\t\"" << representative(myRepeatIterator) << '\"' << ::std::endl;

        ++myRepeatIterator;
    }

    return 0;
}