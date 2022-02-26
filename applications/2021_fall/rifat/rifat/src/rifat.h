#ifndef RIFAT_H
#define RIFAT_H

#include <interfaces/iplugin.h>

class rifat : public KDevelop::IPlugin
{
    Q_OBJECT

public:
    // KPluginFactory-based plugin wants constructor with this signature
    rifat(QObject* parent, const QVariantList& args);
};

#endif // RIFAT_H
