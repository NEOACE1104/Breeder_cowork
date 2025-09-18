import { useTranslation } from 'react-i18next';
import { Card, CardHeader, CardTitle, CardDescription, CardContent } from '@/components/ui/card';

const ReportsPage = () => {
  const { t } = useTranslation();

  return (
    <div className="space-y-6">
      <h1 className="text-3xl font-bold">{t('reports.title')}</h1>
      <Card>
        <CardHeader>
          <CardTitle>{t('reports.title')}</CardTitle>
          <CardDescription>{t('reports.description')}</CardDescription>
        </CardHeader>
        <CardContent>
          <p className="text-muted-foreground">
            Feature under construction. This section will list all your past analysis runs.
          </p>
           <ul className="list-disc list-inside text-muted-foreground mt-2">
            <li>See a history of QC, GP, GWAS, and Mating runs.</li>
            <li>View the parameters used for each run.</li>
            <li>Download previously generated results and HTML reports.</li>
          </ul>
        </CardContent>
      </Card>
    </div>
  );
};

export default ReportsPage;
