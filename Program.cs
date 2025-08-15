namespace Tsumiki
{
    internal class Program
    {
        static void Main(string[] args)
        {
            if (args.Length == 0)
            {
                Console.Error.WriteLine(Consts.DETAILS_TEXT);
            }
            Console.WriteLine("開発中！");
        }
    }
}
