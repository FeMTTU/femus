#include "rifat.h"

#include <debug.h>

#include <KPluginFactory>

K_PLUGIN_FACTORY_WITH_JSON(rifatFactory, "rifat.json", registerPlugin<rifat>(); )

rifat::rifat(QObject *parent, const QVariantList& args)
    : KDevelop::IPlugin(QStringLiteral("rifat"), parent)
{
    Q_UNUSED(args);

    qCDebug(PLUGIN_RIFAT) << "Hello world, my plugin is loaded!";
}

// needed for QObject class created from K_PLUGIN_FACTORY_WITH_JSON
#include "rifat.moc"
