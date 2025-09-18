import { useTranslation } from 'react-i18next';
import { Link } from 'react-router-dom';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { TestTube, Dna, Combine } from 'lucide-react';

const DashboardPage = () => {
  const { t } = useTranslation();

  const features = [
    {
      title: t('nav.data_qc'),
      description: t('data_qc.title'),
      link: '/data-qc',
      icon: <TestTube className="h-8 w-8" />,
    },
    {
      title: t('nav.gp'),
      description: t('gp.description'),
      link: '/gp',
      icon: <Dna className="h-8 w-8" />,
    },
    {
      title: t('nav.gwas'),
      description: t('gwas.description'),
      link: '/gwas',
      icon: <Dna className="h-8 w-8" />,
    },
    {
      title: t('nav.mating'),
      description: t('mating.description'),
      link: '/mating',
      icon: <Combine className="h-8 w-8" />,
    },
  ];

  return (
    <div className="space-y-6">
      <Card className="bg-green-50 dark:bg-green-900/20 border-green-200 dark:border-green-800">
        <CardHeader>
          <CardTitle className="text-3xl text-green-800 dark:text-green-200">{t('dashboard.title')}</CardTitle>
          <CardDescription>{t('dashboard.description')}</CardDescription>
        </CardHeader>
      </Card>

      <div>
        <h2 className="text-2xl font-semibold mb-4">{t('dashboard.quick_links')}</h2>
        <div className="grid gap-6 md:grid-cols-2 lg:grid-cols-4">
          {features.map((feature) => (
            <Card key={feature.title}>
              <CardHeader className="flex flex-row items-center justify-between pb-2">
                <CardTitle className="text-lg font-medium">{feature.title}</CardTitle>
                {feature.icon}
              </CardHeader>
              <CardContent>
                <p className="text-sm text-muted-foreground mb-4">{feature.description}</p>
                <Button asChild>
                  <Link to={feature.link}>Go to {feature.title}</Link>
                </Button>
              </CardContent>
            </Card>
          ))}
        </div>
      </div>
    </div>
  );
};

export default DashboardPage;
