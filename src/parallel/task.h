#pragma once

#include "serialization.h"
#include "ioser.h"

template <class ContextT>
class Task
{
public:
    typedef ContextT ContextType; // used in comper.h
    ContextT context;

    friend ibinstream &operator<<(ibinstream &m, const Task &t)
    {

        m << t.context;
        return m;
    }

    friend obinstream &operator>>(obinstream &m, Task &t)
    {
        m >> t.context;
        return m;
    }

    friend ifbinstream &operator<<(ifbinstream &m, const Task &t)
    {
        m << t.context;
        return m;
    }

    friend ofbinstream &operator>>(ofbinstream &m, Task &t)
    {
        m >> t.context;
        return m;
    }
};